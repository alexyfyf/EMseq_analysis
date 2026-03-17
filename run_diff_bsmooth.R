# Load libraries
library(bsseq)
library(BiocParallel)

# Parse command line arguments
# Usage: Rscript run_diff_bsmooth.R FILE_1 FILE_2 NAME_1 NAME_2 OUT_DIR
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript run_diff_bsmooth.R FILE_1 FILE_2 NAME_1 NAME_2 OUT_DIR")
}

file_1 <- args[1]
file_2 <- args[2]
name_1 <- args[3]
name_2 <- args[4]
output_dir <- args[5]

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 1. Read Bismark files into a BSseq object
files <- c(file_1, file_2)
bs_data <- read.bismark(files, 
                        colData = DataFrame(row.names = c(name_1, name_2)),
                        rmZeroCov = TRUE, 
                        strandCollapse = TRUE)

# Filter for minimum coverage of 10 in both samples
keep <- rowSums(getCoverage(bs_data) >= 10) == 2
bs_data <- bs_data[keep, ]

# 2. Perform Smoothing
# Use 8 cores for parallelization
bs_smoothed <- BSmooth(bs_data, BPPARAM = MulticoreParam(workers = 8))

# 3. Extract smoothed values
# Since we have no replicates, we can't use BSmooth.tstat()
# Instead, we calculate the difference between the two columns
meth_control <- getMeth(bs_smoothed, type = "smooth")[, 1] # name_1
meth_treat <- getMeth(bs_smoothed, type = "smooth")[, 2] # name_2
meth_diff <- meth_treat - meth_control

# 4. Identify Regions of Interest (Heuristic approach)
# Find regions where the difference is > 0.1 (10%)
significant_indices <- which(abs(meth_diff) >= 0.1)

# Extract and save the differentially methylated cytosines/loci (DML)
dml_results <- data.frame(
  chr = as.character(seqnames(bs_smoothed))[significant_indices],
  pos = start(bs_smoothed)[significant_indices],
  meth_control = meth_control[significant_indices],
  meth_treat = meth_treat[significant_indices],
  meth_diff = meth_diff[significant_indices]
)
# Update column names to reflect sample names
colnames(dml_results)[3] <- paste0("meth_", name_1)
colnames(dml_results)[4] <- paste0("meth_", name_2)

# Custom function to group DMLs into DMRs 
find_dmrs <- function(dmls, name1, name2, max_gap = 300, min_cpgs = 3) {
  if (nrow(dmls) == 0) return(data.frame())
  
  # Ensure sorted by chr and pos
  dmls <- dmls[order(dmls$chr, dmls$pos), ]
  
  # Determine direction (hyper or hypo in Treatment vs Control)
  dmls$direction <- ifelse(dmls$meth_diff > 0, "hyper", "hypo")
  
  # Calculate distance to previous DML
  dmls$dist <- c(NA, diff(dmls$pos))
  
  # A new cluster starts if: 
  # 1. New chromosome
  # 2. Distance > max_gap
  # 3. Direction changes
  new_cluster <- c(TRUE, dmls$chr[-1] != dmls$chr[-nrow(dmls)] | 
                         dmls$dist[-1] > max_gap |
                         dmls$direction[-1] != dmls$direction[-nrow(dmls)])
  
  # Assign cluster IDs
  dmls$cluster_id <- cumsum(new_cluster)
  
  # Aggregate by cluster
  dmrs <- aggregate(pos ~ cluster_id, data = dmls, FUN = function(x) c(start = min(x), end = max(x)))
  
  # Columns to aggregate
  col1 <- paste0("meth_", name1)
  col2 <- paste0("meth_", name2)
  
  dmrs_info <- aggregate(as.formula(paste("cbind(meth_diff,", col1, ",", col2, ") ~ cluster_id")), data = dmls, FUN = mean)
  chr_dir <- aggregate(cbind(chr, direction) ~ cluster_id, data = dmls, FUN = function(x) x[1])
  counts <- aggregate(pos ~ cluster_id, data = dmls, FUN = length)
  colnames(counts)[2] <- "n_cpgs"
  
  # Combine results
  res <- data.frame(
    chr = chr_dir$chr,
    start = dmrs$pos[, "start"],
    end = dmrs$pos[, "end"],
    n_cpgs = counts$n_cpgs,
    mean_meth_diff = dmrs_info$meth_diff,
    mean_meth_control = dmrs_info[[col1]],
    mean_meth_treat = dmrs_info[[col2]],
    direction = chr_dir$direction
  )
  colnames(res)[6] <- paste0("mean_meth_", name1)
  colnames(res)[7] <- paste0("mean_meth_", name2)
  
  # Filter by minimum CpGs per DMR
  res <- res[res$n_cpgs >= min_cpgs, ]
  return(res)
}

dmr_results <- find_dmrs(dml_results, name_1, name_2, max_gap = 300, min_cpgs = 3)

# Save results
write.csv(dml_results, file=file.path(output_dir, "bsmooth_sig_diff_meth_10p_loci.csv"), row.names=FALSE)
write.csv(dmr_results, file=file.path(output_dir, "bsmooth_sig_diff_meth_10p_regions.csv"), row.names=FALSE)

print("BSmooth analysis complete.")