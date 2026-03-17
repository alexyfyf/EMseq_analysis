library(methylKit)
library(annotatr)
library(GenomicRanges)

output_dir = "/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq/data"

print("Loading data...")
data <- readRDS(file.path(output_dir, "methylKit_unite.rds"))
meth <- data$meth

# Annotate with CpG islands, shores, shelves (and intergenic)
print("Annotating with annotatr...")
# Get annotations for hg38. cpg_inter is already included in hg38_cpgs
cpgs_info <- build_annotations(genome = "hg38", annotations = "hg38_cpgs")

# Convert unite object to GRanges to overlap with annotations
# We will use the coverage data (total number of reads per base)
print("Calculating read coverage over annotations...")
meth_gr <- as(meth, "GRanges")
cvg_A <- mcols(meth_gr)$coverage1 # Assuming coverage1 is CB_A_1 and coverage2 is CB_V_1
cvg_V <- mcols(meth_gr)$coverage2 
mcols(meth_gr)$cvg_A <- cvg_A
mcols(meth_gr)$cvg_V <- cvg_V

# Intersect methylation data with annotations
overlaps <- findOverlaps(meth_gr, cpgs_info)

# Create a data frame summarizing the overlaps for coverage
meth_anno <- data.frame(
  cvg_A = mcols(meth_gr[queryHits(overlaps)])$cvg_A,
  cvg_V = mcols(meth_gr[queryHits(overlaps)])$cvg_V,
  annot_type = as.character(cpgs_info[subjectHits(overlaps)]$type)
)

# Sum of all read bases in each annotation type
sum_cvg_A <- aggregate(cvg_A ~ annot_type, data = meth_anno, FUN = sum, na.rm=TRUE)
sum_cvg_V <- aggregate(cvg_V ~ annot_type, data = meth_anno, FUN = sum, na.rm=TRUE)

# To get total length of each annotation type, we merge overlapping regions of the same type
annot_lengths <- sapply(split(cpgs_info, cpgs_info$type), function(gr) {
  sum(width(reduce(gr)))
})
annot_lengths_df <- data.frame(annot_type = names(annot_lengths), total_bp = as.numeric(annot_lengths))

# Combine sums and lengths
summary_cpgs <- merge(sum_cvg_A, sum_cvg_V, by="annot_type", all=TRUE)
summary_cpgs <- merge(summary_cpgs, annot_lengths_df, by="annot_type", all.x=TRUE)

# Calculate coverage ratios (total reads / total bp)
summary_cpgs$cov_ratio_A <- summary_cpgs$cvg_A / summary_cpgs$total_bp
summary_cpgs$cov_ratio_V <- summary_cpgs$cvg_V / summary_cpgs$total_bp

summary_cpgs = summary_cpgs[,c("annot_type", "total_bp", "cvg_A", "cvg_V", "cov_ratio_A", "cov_ratio_V")]

# Save the results
print("Saving results...")
write.csv(summary_cpgs, file=file.path(output_dir, "methylKit_CpG_regional_summary.csv"), row.names=FALSE)

print("Analysis complete.")
