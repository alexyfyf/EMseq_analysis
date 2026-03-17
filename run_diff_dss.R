# Load libraries
library(DSS)
library(bsseq)

# Parse command line arguments
# Usage: Rscript run_diff_dss.R FILE_1 FILE_2 NAME_1 NAME_2 OUT_DIR
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript run_diff_dss.R FILE_1 FILE_2 NAME_1 NAME_2 OUT_DIR")
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

# 1. Load Bismark coverage files
# DSS requires a data frame with columns: chr, pos, N, X
# N = total reads, X = methylated reads
sample1 <- read.table(file_1, header=FALSE)
sample2 <- read.table(file_2, header=FALSE)

colnames(sample1) <- c("chr", "pos", "end", "perc", "X", "U")
colnames(sample2) <- c("chr", "pos", "end", "perc", "X", "U")

# Create N (total count)
sample1$N <- sample1$X + sample1$U
sample2$N <- sample2$X + sample2$U

# Filter for minimum coverage of 10
sample1 <- subset(sample1, N >= 10)
sample2 <- subset(sample2, N >= 10)

# 2. Format for DSS (only need chr, pos, N, X)
dat1 <- sample1[, c("chr", "pos", "N", "X")]
dat2 <- sample2[, c("chr", "pos", "N", "X")]

# 3. Create a BSseq object
BSobj <- makeBSseqData(list(dat1, dat2), c(name_1, name_2))

# 4. Perform DML (Differential Methylation Loci) test
# Smoothing is CRITICAL when you have no replicates
dmlTest <- DMLtest(BSobj, group1=name_1, group2=name_2, smoothing=TRUE)

# Call explicit differentially methylated loci (DML)
dmls <- callDML(dmlTest, delta=0.1, p.threshold=0.01)

# 5. Call DMRs (Regions)
dmrs <- callDMR(dmlTest, delta=0.1, p.threshold=0.01)

# Save output
write.csv(dmls, file=file.path(output_dir, "dss_sig_diff_meth_10p_loci.csv"), row.names=FALSE)
write.csv(dmrs, file=file.path(output_dir, "dss_sig_diff_meth_10p_regions.csv"), row.names=FALSE)

print("DSS analysis complete.")