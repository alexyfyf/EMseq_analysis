library(tidyverse)
library(UpSetR)
library(ChIPpeakAnno)
library(GenomicRanges)
library(ggalluvial)
library(plyranges)

# Source utility functions
source("/home/users/allstaff/yan.a/yan.a/software/utils/join_overlap_full.R")
source("/home/users/allstaff/yan.a/yan.a/software/utils/plot_sankey.R")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript compare_cell_line_dss.R CELL_LINE OUT_DIR")
}

cell_line <- args[1]
out_dir <- args[2]

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Base directory for results
base_output_dir <- "/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq/output"

# Define available comparisons
comparison_configs <- list(
  "A_vs_G" = list(dir = file.path(base_output_dir, paste0(cell_line, "_A_vs_G")), target = "A_vs_G"),
  "A_vs_V" = list(dir = file.path(base_output_dir, paste0(cell_line, "_A_vs_V")), target = "A_vs_V"),
  "G_vs_V" = list(dir = file.path(base_output_dir, paste0(cell_line, "_G_vs_V")), target = "G_vs_V")
)

# Define comparison subsets to plot
comparison_subsets <- list(
  "3way" = c("A_vs_G", "A_vs_V", "G_vs_V"),
  "2way" = c("A_vs_V", "G_vs_V")
)

# Function to read and normalize comparison results
load_and_normalize <- function(comp_name, config, type = "dml") {
  dir_path <- config$dir
  file_name <- if (type == "dml") "dss_sig_diff_meth_10p_loci.csv" else "dss_sig_diff_meth_10p_regions.csv"
  file_path <- file.path(dir_path, file_name)
  
  if (!file.exists(file_path)) {
    message(paste("File missing:", file_path))
    return(NULL)
  }
  
  df <- read.csv(file_path)
  if (nrow(df) == 0) return(NULL)
  
  # Standardize direction (Keep as delivered in folders)
  df$norm_diff <- df$diff
  df$target_comp <- config$target
  
  # Convert to GRanges
  if (type == "dml") {
    gr <- GRanges(seqnames = df$chr, ranges = IRanges(start = df$pos, end = df$pos))
  } else {
    gr <- GRanges(seqnames = df$chr, ranges = IRanges(start = df$start, end = df$end))
  }
  mcols(gr) <- df[, c("norm_diff", "target_comp")]
  
  return(gr)
}

# General processing function
process_comparisons <- function(type, subset_name, comp_names) {
  message(sprintf("Processing %s for %s subset...", type, subset_name))
  
  results_list <- list()
  for (name in comp_names) {
    gr <- load_and_normalize(name, comparison_configs[[name]], type = type)
    if (!is.null(gr)) {
      results_list[[name]] <- gr
    }
  }
  
  if (length(results_list) == 0) {
    message("No data found for this subset.")
    return(NULL)
  }
  
  pdf_file <- file.path(out_dir, paste0(cell_line, "_", type, "_", subset_name, ".pdf"))
  pdf(pdf_file, width = 12, height = 8)
  
  # 1. UpSet or Venn Plot
  if (type == "dml") {
    hyper_sets <- lapply(results_list, function(x) paste(seqnames(x), start(x), sep="_")[mcols(x)$norm_diff > 0])
    hypo_sets <- lapply(results_list, function(x) paste(seqnames(x), start(x), sep="_")[mcols(x)$norm_diff < 0])
    
    if (length(unlist(hyper_sets)) > 0) {
      print(upset(fromList(hyper_sets), order.by = "freq"))
      grid::grid.text(paste(cell_line, "DML Hyper-methylated", subset_name), x = 0.65, y = 0.95, gp = grid::gpar(fontsize = 15))
    }
    
    if (length(unlist(hypo_sets)) > 0) {
      print(upset(fromList(hypo_sets), order.by = "freq"))
      grid::grid.text(paste(cell_line, "DML Hypo-methylated", subset_name), x = 0.65, y = 0.95, gp = grid::gpar(fontsize = 15))
    }
  } else {
    # DMRs use ChIPpeakAnno
    hyper_gr <- lapply(results_list, function(x) x[mcols(x)$norm_diff > 0])
    hypo_gr <- lapply(results_list, function(x) x[mcols(x)$norm_diff < 0])
    
    if (length(hyper_gr) >= 2) {
      makeVennDiagram(hyper_gr, NameOfPeaks = names(hyper_gr), main = paste(cell_line, "DMR Hyper Overlap", subset_name))
    }
    
    if (length(hypo_gr) >= 2) {
      makeVennDiagram(hypo_gr, NameOfPeaks = names(hypo_gr), main = paste(cell_line, "DMR Hypo Overlap", subset_name))
    }
  }
  
  # 2. Sankey Plot
  if (length(results_list) >= 2) {
    message("Generating Sankey plot...")
    plot_sankey(results_list, column = "norm_diff")
    grid::grid.text(paste(cell_line, type, "Sankey", subset_name), x = 0.5, y = 0.95, gp = grid::gpar(fontsize = 15))
  }
  
  dev.off()
  message(paste("Saved", pdf_file))
}

# Run for all combinations
for (subset_name in names(comparison_subsets)) {
  comp_subset <- comparison_subsets[[subset_name]]
  process_comparisons("dml", subset_name, comp_subset)
  process_comparisons("dmr", subset_name, comp_subset)
}
