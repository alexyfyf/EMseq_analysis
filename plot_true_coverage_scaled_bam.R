library(genomation)
library(annotatr)
library(GenomicRanges)

# 1. Define paths
bam_file <- "/home/users/allstaff/yan.a/davidson_longread/yan.a/XZ/emseq/results_Pt_A/bismark/deduplicated/Pt_A.deduplicated.sorted.bam"
output_pdf <- "/home/users/allstaff/yan.a/davidson_longread/yan.a/XZ/emseq/data/Pt_A_true_coverage_SCALED.pdf"

# 2. Get CpG annotations for hg38
message("Fetching hg38 CpG annotations...")
annots <- c('hg38_cpg_islands', 'hg38_cpg_shores', 'hg38_cpg_shelves')
annotations <- build_annotations(genome = 'hg38', annotations = annots)
annotations <- keepStandardChromosomes(annotations, pruning.mode="coarse")

# 3. Split annotations by type
split_annots <- split(annotations, annotations$type)

# 4. Calculate Score Matrices (Scaling each region to 100 bins)
message("Calculating scaled coverage from BAM...")
score_list <- lapply(names(split_annots), function(name) {
  message("  Processing ", name, "...")
  regions <- split_annots[[name]]
  
  # Filter for regions with width >= 100 to allow division into 100 bins
  regions <- regions[width(regions) >= 100]
  
  # Use ScoreMatrixBin on the full regions to scale them to 100 bins
  sm <- tryCatch({
    ScoreMatrixBin(target = bam_file, windows = regions, bin.num = 100, type = "bam")
  }, error = function(e) {
    message("    Error: ", e$message)
    return(NULL)
  })
  return(sm)
})

# Filter out NULLs and assign names
score_list <- score_list[!vapply(score_list, is.null, logical(1))]
names(score_list) <- names(split_annots)[names(split_annots) %in% names(score_list)]

# 5. Plotting
if (length(score_list) > 0) {
  message("Plotting...")
  pdf(output_pdf, width = 10, height = 7)
  
  labels <- c("CpG Islands", "CpG Shelves", "CpG Shores")
  cols <- c("#1b9e77", "#d95f02", "#7570b3")
  
  plotMeta(ScoreMatrixList(score_list), 
           overlay = TRUE, 
           centralTend = "mean",
           profile.names = labels,
           main = "True Coverage Profile (Scaled 0-100% of Region Width)",
           xlab = "Relative Position (Start to End of Region)",
           ylab = "Mean Read Depth (Coverage)",
           line.col = cols)
  
  legend("topright", legend = labels, col = cols, lty = 1, lwd = 2, bty = "n")
  
  dev.off()
  message("Done! Plot saved to: ", output_pdf)
}
