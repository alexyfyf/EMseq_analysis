#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(annotatr)
  library(genomation)
})

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript plot_cpg_region_coverage.R --input <sample.bismark.cov.gz>",
      "    [--genome hg38] [--out-prefix sample]",
      "    [--chunk-size 200000] [--sample-frac 0.02] [--bins 100]",
      "",
      "Required packages:",
      "  data.table, GenomicRanges, IRanges, annotatr, genomation",
      "",
      "Notes:",
      "  - Coverage is methylated + unmethylated reads from the Bismark .cov.gz file.",
      "  - CpG classes come from annotatr built-in CpG annotations.",
      "  - The genomation meta-profile uses a random sample of loci for tractable plotting.",
      "  - Genomation plots are centered on each region and extended +/- 1 kb.",
      "  - The script writes per-covered-CpG coverage and per-region mean coverage tables.",
      sep = "\n"
    )
  )
}

parse_args <- function(argv) {
  defaults <- list(
    input = NULL,
    genome = "hg38",
    out_prefix = NULL,
    chunk_size = 200000L,
    sample_frac = 0.02,
    bins = 100L
  )

  if (length(argv) == 0L || any(argv %in% c("-h", "--help"))) {
    usage()
    quit(save = "no", status = 0L)
  }

  i <- 1L
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    if (i == length(argv)) {
      stop("Missing value for ", key, call. = FALSE)
    }
    value <- argv[[i + 1L]]
    name <- sub("^--", "", key)
    name <- chartr("-", "_", name)
    if (!name %in% names(defaults)) {
      stop("Unknown argument: ", key, call. = FALSE)
    }
    defaults[[name]] <- value
    i <- i + 2L
  }

  defaults$chunk_size <- as.integer(defaults$chunk_size)
  defaults$sample_frac <- as.numeric(defaults$sample_frac)
  defaults$bins <- as.integer(defaults$bins)

  if (is.null(defaults$input)) {
    stop("--input is required", call. = FALSE)
  }
  if (is.null(defaults$out_prefix)) {
    stem <- basename(defaults$input)
    stem <- sub("\\.gz$", "", stem)
    stem <- sub("\\.bismark\\.cov$", "", stem)
    defaults$out_prefix <- stem
  }
  if (defaults$sample_frac <= 0 || defaults$sample_frac > 1) {
    stop("--sample-frac must be in the interval (0, 1]", call. = FALSE)
  }
  if (defaults$chunk_size < 1L) {
    stop("--chunk-size must be >= 1", call. = FALSE)
  }
  if (defaults$bins < 5L) {
    stop("--bins must be >= 5", call. = FALSE)
  }

  defaults
}

normalize_cpg_class <- function(x) {
  x <- tolower(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("island", x)] <- "island"
  out[grepl("shore", x)] <- "shore"
  out[grepl("shelf", x)] <- "shelf"
  out[grepl("inter|open", x)] <- "open_sea"
  out
}

pretty_cpg_class <- function(x) {
  labels <- c(
    island = "CpG island",
    shore = "Shore",
    shelf = "Shelf",
    open_sea = "Open sea"
  )
  unname(labels[x])
}

read_cov_chunk <- function(con, n) {
  lines <- readLines(con, n = n)
  if (length(lines) == 0L) {
    return(NULL)
  }

  dt <- fread(
    text = lines,
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "methylation_pct", "methylated", "unmethylated"),
    showProgress = FALSE
  )

  dt[, coverage := methylated + unmethylated]
  dt
}

extract_annotation_labels <- function(gr) {
  cols <- colnames(S4Vectors::mcols(gr))
  candidates <- c("type", "annot.type", "annotation", "id")

  for (col in candidates[candidates %in% cols]) {
    vals <- S4Vectors::mcols(gr)[[col]]

    if (length(vals) == length(gr)) {
      return(as.character(vals))
    }

    if (methods::is(vals, "CompressedCharacterList") || methods::is(vals, "CharacterList")) {
      flat <- vapply(as.list(vals), function(x) {
        if (length(x) == 0L) {
          NA_character_
        } else {
          paste(x, collapse = ";")
        }
      }, character(1))
      if (length(flat) == length(gr)) {
        return(flat)
      }
    }
  }

  nm <- names(gr)
  if (!is.null(nm) && length(nm) == length(gr)) {
    return(nm)
  }

  stop("Could not extract one annotation label per range from annotatr output.", call. = FALSE)
}

build_cpg_annotations <- function(genome) {
  annotation_map <- c(
    island = paste0(genome, "_cpg_islands"),
    shore = paste0(genome, "_cpg_shores"),
    shelf = paste0(genome, "_cpg_shelves"),
    open_sea = paste0(genome, "_cpg_inter")
  )

  annos <- annotatr::builtin_annotations()
  missing <- annotation_map[!annotation_map %in% annos]
  if (length(missing) > 0L) {
    stop("Missing annotatr CpG annotations: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  gr_list <- lapply(names(annotation_map), function(region_name) {
    gr <- annotatr::build_annotations(genome = genome, annotations = annotation_map[[region_name]])
    S4Vectors::mcols(gr)$cpg_class <- region_name
    gr
  })
  names(gr_list) <- names(annotation_map)

  annot_gr <- do.call(c, unname(gr_list))
  S4Vectors::mcols(annot_gr)$region_index <- seq_along(annot_gr)
  annot_gr
}

make_score_matrix <- function(target_gr, windows_gr, bins) {
  if (length(target_gr) == 0L || length(windows_gr) == 0L) {
    return(NULL)
  }

  genomation::ScoreMatrixBin(
    target = target_gr,
    windows = windows_gr,
    bin.num = bins,
    bin.op = "mean",
    weight.col = "score",
    strand.aware = FALSE
  )
}

make_centered_windows <- function(gr, flank = 1000L) {
  if (length(gr) == 0L) {
    return(gr)
  }

  centers <- start(gr) + (width(gr) - 1L) %/% 2L
  GRanges(
    seqnames = seqnames(gr),
    ranges = IRanges(start = pmax(1L, centers - flank), end = centers + flank),
    strand = strand(gr),
    S4Vectors::mcols(gr)
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
set.seed(1L)

message("Building annotatr CpG annotations for ", args$genome, " ...")
annot_gr <- build_cpg_annotations(args$genome)

region_levels <- c("island", "shore", "shelf", "open_sea")
region_labels <- pretty_cpg_class(region_levels)
names(region_labels) <- region_levels

annot_by_class <- setNames(
  lapply(region_levels, function(region_name) {
    annot_gr[S4Vectors::mcols(annot_gr)$cpg_class == region_name]
  }),
  region_levels
)
region_count_dt <- data.table(
  region = region_levels,
  annotated_regions = vapply(annot_by_class, length, integer(1)),
  width_gt1_regions = vapply(annot_by_class, function(gr) sum(width(gr) > 1L), integer(1))
)

stats_dt <- data.table(
  region = region_levels,
  sites = 0L,
  total_coverage = 0
)

sample_cov <- setNames(vector("list", length(region_levels)), region_levels)
full_cov <- setNames(vector("list", length(region_levels)), region_levels)
sample_gr <- GRanges()
meta_gr <- GRanges()
region_sum_coverage <- numeric(length(annot_gr))
region_cpg_count <- integer(length(annot_gr))

per_cpg_path <- paste0(args$out_prefix, "_per_cpg_coverage.tsv.gz")
fwrite(
  data.table(
    chr = character(),
    start = integer(),
    end = integer(),
    coverage = integer(),
    cpg_category = character()
  ),
  per_cpg_path,
  sep = "\t"
)

message("Reading coverage in chunks from ", args$input, " ...")
con <- gzfile(args$input, open = "rt")
on.exit(close(con), add = TRUE)

chunk_id <- 0L
repeat {
  chunk <- read_cov_chunk(con, args$chunk_size)
  if (is.null(chunk)) {
    break
  }
  chunk_id <- chunk_id + 1L

  cov_gr <- GRanges(
    seqnames = chunk$chr,
    ranges = IRanges(start = chunk$start, end = chunk$end),
    score = chunk$coverage
  )

  hits <- findOverlaps(cov_gr, annot_gr, ignore.strand = TRUE)
  if (length(hits) == 0L) {
    next
  }

  query_idx <- queryHits(hits)
  subject_idx <- subjectHits(hits)
  matched_regions <- S4Vectors::mcols(annot_gr)$cpg_class[subject_idx]
  matched_cov <- S4Vectors::mcols(cov_gr)$score[query_idx]
  matched_region_index <- S4Vectors::mcols(annot_gr)$region_index[subject_idx]

  region_accum_dt <- data.table(region_index = matched_region_index, coverage = matched_cov)[
    ,
    .(sum_coverage = sum(coverage), cpg_count = .N),
    by = region_index
  ]
  region_sum_coverage[region_accum_dt$region_index] <- region_sum_coverage[region_accum_dt$region_index] + region_accum_dt$sum_coverage
  region_cpg_count[region_accum_dt$region_index] <- region_cpg_count[region_accum_dt$region_index] + region_accum_dt$cpg_count

  chunk_stats <- data.table(region = matched_regions, coverage = matched_cov)[
    ,
    .(sites = .N, total_coverage = sum(coverage)),
    by = region
  ]
  stats_dt[chunk_stats, on = "region", `:=`(
    sites = sites + i.sites,
    total_coverage = total_coverage + i.total_coverage
  )]

  full_split_vals <- split(matched_cov, matched_regions)
  for (nm in names(full_split_vals)) {
    full_cov[[nm]] <- c(full_cov[[nm]], full_split_vals[[nm]])
  }

  per_cpg_dt <- unique(data.table(
    chr = as.character(seqnames(cov_gr)[query_idx]),
    start = start(cov_gr)[query_idx],
    end = end(cov_gr)[query_idx],
    coverage = matched_cov,
    cpg_category = region_labels[matched_regions]
  ))
  fwrite(per_cpg_dt, per_cpg_path, sep = "\t", append = TRUE, col.names = FALSE)

  meta_gr <- c(meta_gr, cov_gr[unique(query_idx)])

  sampled_rows <- which(stats::runif(nrow(chunk)) <= args$sample_frac)
  if (length(sampled_rows) > 0L) {
    sampled_gr <- cov_gr[sampled_rows]
    sample_gr <- c(sample_gr, sampled_gr)

    sampled_hits <- findOverlaps(sampled_gr, annot_gr, ignore.strand = TRUE)
    if (length(sampled_hits) > 0L) {
      sampled_query <- queryHits(sampled_hits)
      sampled_subject <- subjectHits(sampled_hits)
      sampled_regions <- S4Vectors::mcols(annot_gr)$cpg_class[sampled_subject]
      sampled_cov_values <- S4Vectors::mcols(sampled_gr)$score[sampled_query]
      split_vals <- split(sampled_cov_values, sampled_regions)
      for (nm in names(split_vals)) {
        sample_cov[[nm]] <- c(sample_cov[[nm]], split_vals[[nm]])
      }
    }
  }

  if (chunk_id %% 25L == 0L) {
    message("Processed ", format(chunk_id, big.mark = ","), " chunks")
  }
}

stats_dt[, mean_coverage := fifelse(sites > 0L, total_coverage / sites, NA_real_)]
stats_dt[, region_label := region_labels[region]]

summary_path <- paste0(args$out_prefix, "_cpg_region_coverage_summary.tsv")
fwrite(
  stats_dt[, .(region = region_label, sites, total_coverage, mean_coverage)],
  summary_path,
  sep = "\t"
)

region_mean_dt <- data.table(
  chr = as.character(seqnames(annot_gr)),
  start = start(annot_gr),
  end = end(annot_gr),
  region = region_labels[S4Vectors::mcols(annot_gr)$cpg_class],
  covered_cpg_sites = region_cpg_count,
  total_coverage = region_sum_coverage,
  mean_coverage = fifelse(region_cpg_count > 0L, region_sum_coverage / region_cpg_count, NA_real_)
)
region_mean_path <- paste0(args$out_prefix, "_per_region_mean_coverage.tsv")
fwrite(region_mean_dt, region_mean_path, sep = "\t")

boxplot_df <- rbindlist(
  lapply(region_levels, function(region_name) {
    vals <- sample_cov[[region_name]]
    if (length(vals) == 0L) {
      return(data.table(region = region_labels[[region_name]], coverage = NA_real_))
    }
    data.table(region = region_labels[[region_name]], coverage = vals)
  }),
  use.names = TRUE
)

if (nrow(boxplot_df) > 0L) {
  png(
    filename = paste0(args$out_prefix, "_coverage_boxplot_sampled.png"),
    width = 1400,
    height = 900,
    res = 150
  )
  boxplot(
    coverage ~ region,
    data = boxplot_df,
    log = "y",
    col = c("#1b9e77", "#d95f02", "#7570b3", "#666666"),
    outline = FALSE,
    las = 2,
    ylab = "Coverage per CpG",
    xlab = "",
    main = "Coverage by CpG context (sampled, log-scaled y-axis)"
  )
  dev.off()
}

boxplot_full_df <- rbindlist(
  lapply(region_levels, function(region_name) {
    vals <- full_cov[[region_name]]
    if (length(vals) == 0L) {
      return(data.table(region = region_labels[[region_name]], coverage = NA_real_))
    }
    data.table(region = region_labels[[region_name]], coverage = vals)
  }),
  use.names = TRUE
)

if (nrow(boxplot_full_df) > 0L) {
  png(
    filename = paste0(args$out_prefix, "_coverage_boxplot_full.png"),
    width = 1400,
    height = 900,
    res = 150
  )
  boxplot(
    coverage ~ region,
    data = boxplot_full_df,
    log = "y",
    col = c("#1b9e77", "#d95f02", "#7570b3", "#666666"),
    outline = FALSE,
    las = 2,
    ylab = "Coverage per CpG",
    xlab = "",
    main = "Coverage by CpG context (all covered CpGs, log-scaled y-axis)"
  )
  dev.off()
}

message("Building genomation score matrices on centered +/-1 kb windows ...")
score_mats_sampled <- lapply(region_levels, function(region_name) {
  windows <- make_centered_windows(annot_by_class[[region_name]], flank = 1000L)
  make_score_matrix(sample_gr, windows, bins = args$bins)
})
names(score_mats_sampled) <- region_levels
score_mats_sampled <- score_mats_sampled[!vapply(score_mats_sampled, is.null, logical(1))]

score_mats_full <- lapply(region_levels, function(region_name) {
  windows <- make_centered_windows(annot_by_class[[region_name]], flank = 1000L)
  make_score_matrix(meta_gr, windows, bins = args$bins)
})
names(score_mats_full) <- region_levels
score_mats_full <- score_mats_full[!vapply(score_mats_full, is.null, logical(1))]

plot_usage_dt <- data.table(
  region = region_levels,
  plotted_regions_sampled = vapply(region_levels, function(region_name) {
    sm <- score_mats_sampled[[region_name]]
    if (is.null(sm)) {
      return(0L)
    }
    nrow(sm)
  }, integer(1)),
  plotted_regions_full = vapply(region_levels, function(region_name) {
    sm <- score_mats_full[[region_name]]
    if (is.null(sm)) {
      return(0L)
    }
    nrow(sm)
  }, integer(1)),
  sampled_loci_for_boxplot = vapply(region_levels, function(region_name) {
    length(sample_cov[[region_name]])
  }, integer(1)),
  meta_loci_sampled = length(sample_gr),
  meta_loci_full = length(meta_gr)
)

report_dt <- merge(region_count_dt, plot_usage_dt, by = "region", all = TRUE)
report_dt[, region_label := region_labels[region]]
report_path <- paste0(args$out_prefix, "_cpg_region_plotting_counts.tsv")
fwrite(
  report_dt[, .(region = region_label, annotated_regions, width_gt1_regions, plotted_regions_sampled, plotted_regions_full, sampled_loci_for_boxplot, meta_loci_sampled, meta_loci_full)],
  report_path,
  sep = "\t"
)

message("CpG category counts:")
for (i in seq_len(nrow(report_dt))) {
  message(
    "  ", report_dt$region_label[[i]],
    ": annotated_regions=", format(report_dt$annotated_regions[[i]], big.mark = ","),
    ", width_gt1_regions=", format(report_dt$width_gt1_regions[[i]], big.mark = ","),
    ", plotted_regions_sampled=", format(report_dt$plotted_regions_sampled[[i]], big.mark = ","),
    ", plotted_regions_full=", format(report_dt$plotted_regions_full[[i]], big.mark = ","),
    ", sampled_loci_for_boxplot=", format(report_dt$sampled_loci_for_boxplot[[i]], big.mark = ","),
    ", meta_loci_sampled=", format(report_dt$meta_loci_sampled[[i]], big.mark = ","),
    ", meta_loci_full=", format(report_dt$meta_loci_full[[i]], big.mark = ",")
  )
}

if (length(score_mats_sampled) > 0L) {
  sml <- methods::new("ScoreMatrixList", score_mats_sampled)
  png(
    filename = paste0(args$out_prefix, "_genomation_meta_profile_sampled.png"),
    width = 1400,
    height = 900,
    res = 150
  )
  plotMeta(
    mat = sml,
    overlay = TRUE,
    profile.names = region_labels[names(score_mats_sampled)],
    centralTend = "mean",
    xlab = "Position relative to region center (-1 kb to +1 kb)",
    ylab = "Average per-CpG coverage at this relative position around the region center",
    line.col = c("#1b9e77", "#d95f02", "#7570b3", "#666666")[seq_along(score_mats_sampled)],
    main = "Coverage around centered CpG-context regions (sampled)"
  )
  dev.off()
}

if (length(score_mats_full) > 0L) {
  sml <- methods::new("ScoreMatrixList", score_mats_full)
  png(
    filename = paste0(args$out_prefix, "_genomation_meta_profile_full.png"),
    width = 1400,
    height = 900,
    res = 150
  )
  plotMeta(
    mat = sml,
    overlay = TRUE,
    profile.names = region_labels[names(score_mats_full)],
    centralTend = "mean",
    xlab = "Position relative to region center (-1 kb to +1 kb)",
    ylab = "Average per-CpG coverage at this relative position around the region center",
    line.col = c("#1b9e77", "#d95f02", "#7570b3", "#666666")[seq_along(score_mats_full)],
    main = "Coverage around centered CpG-context regions (full data)"
  )
  dev.off()
}

message("Wrote:")
message("  ", summary_path)
message("  ", per_cpg_path)
message("  ", region_mean_path)
message("  ", report_path)
if (nrow(boxplot_df) > 0L) {
  message("  ", paste0(args$out_prefix, "_coverage_boxplot_sampled.png"))
}
if (nrow(boxplot_full_df) > 0L) {
  message("  ", paste0(args$out_prefix, "_coverage_boxplot_full.png"))
}
if (length(score_mats_sampled) > 0L) {
  message("  ", paste0(args$out_prefix, "_genomation_meta_profile_sampled.png"))
}
if (length(score_mats_full) > 0L) {
  message("  ", paste0(args$out_prefix, "_genomation_meta_profile_full.png"))
}
