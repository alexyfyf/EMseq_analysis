library(genomation)
library(annotatr)
library(GenomicRanges)

# 1. Define paths
bam_file <- "/home/users/allstaff/yan.a/davidson_longread/yan.a/XZ/emseq/results_Pt_A/bismark/deduplicated/Pt_A.deduplicated.sorted.bam"
output_pdf <- "/home/users/allstaff/yan.a/davidson_longread/yan.a/XZ/emseq/data/Pt_A_true_genomic_coverage.pdf"

# 2. Get CpG annotations for hg38
message("Fetching hg38 CpG annotations...")
annots <- c('hg38_cpg_islands', 'hg38_cpg_shores', 'hg38_cpg_shelves')
annotations <- build_annotations(genome = 'hg38', annotations = annots)

# Filter for standard chromosomes only
annotations <- keepStandardChromosomes(annotations, pruning.mode="coarse")

# 3. Prepare centered windows (+/- 1kb around the center of each region)
message("Preparing windows...")
make_windows <- function(gr) {
  # Remove regions that are too small or near chromosome edges to avoid errors
  gr <- gr[width(gr) > 0]
  centers <- resize(gr, width = 1, fix = "center")
  windows <- promoters(centers, upstream = 1000, downstream = 1000)
  # Trim windows to ensure they stay within chromosome boundaries
  windows <- trim(windows)
  # Keep only windows that are exactly 2000 bp (remove those trimmed at edges)
  windows <- windows[width(windows) == 2000]
  return(windows)
}

# Split annotations by type
split_annots <- split(annotations, annotations$type)

# 4. Calculate Score Matrices from BAM (this counts actual reads)
message("Calculating coverage from BAM (this may take a few minutes)...")
score_list <- lapply(names(split_annots), function(name) {
  message("  Processing ", name, "...")
  windows <- make_windows(split_annots[[name]])
  
  # Try to calculate score matrix, return NULL if it fails
  sm <- tryCatch({
    ScoreMatrixBin(target = bam_file, windows = windows, bin.num = 100, type = "bam")
  }, error = function(e) {
    message("    Error processing ", name, ": ", e$message)
    return(NULL)
  })
  return(sm)
})

# Filter out any NULL matrices
score_list <- score_list[!vapply(score_list, is.null, logical(1))]
names(score_list) <- names(split_annots)[names(split_annots) %in% names(score_list)]

# 5. Plot the result
if (length(score_list) > 0) {
  message("Plotting...")
  pdf(output_pdf, width = 10, height = 7)
  
  # Define nice labels for the legend
  labels <- c("CpG Islands", "CpG Shelves", "CpG Shores")
  cols <- c("#1b9e77", "#d95f02", "#7570b3")
  
  plotMeta(ScoreMatrixList(score_list), 
           overlay = TRUE, 
           centralTend = "mean",
           profile.names = labels,
           main = "True Genomic Coverage (Read Depth) around CpG Regions",
           xlab = "Relative Position (Center +/- 1kb)",
           ylab = "Mean Read Depth (Coverage)",
           line.col = cols)
  
  # Add the legend manually (genomation plotMeta doesn't add one by default when overlay=TRUE)
  legend("topright", legend = labels, col = cols, lty = 1, lwd = 2, bty = "n")
  
  dev.off()
  message("Done! Plot saved to: ", output_pdf)
} else {
  message("No data to plot.")
}
