library(methylKit)

# Parse command line arguments
# Usage: Rscript run_diff_methylkit.R FILE_1 FILE_2 NAME_1 NAME_2 OUT_DIR
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript run_diff_methylkit.R FILE_1 FILE_2 NAME_1 NAME_2 OUT_DIR")
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

# Define file paths
file.list = list(file_1, file_2)

# Read the data matching user's standard parameters
print(paste("Reading methylation data for", name_1, "and", name_2))
myobj = methRead(
  file.list,
  sample.id = list(name_1, name_2),
  assembly = "hg38",
  treatment = c(1, 0),
  context = "CpG",
  pipeline = "bismarkCoverage",
  header = FALSE,
  mincov = 10
)

# Merge samples
print("Merging samples...")
meth = unite(myobj, destrand = FALSE, min.per.group = 1L)

# Calculate differential methylation
print("Calculating differential methylation...")
myDiff = calculateDiffMeth(meth, mc.cores = 8) 

# Get significantly differentially methylated bases
# Q-value < 0.01 and diff methylation > 10%
print("Filtering for significant differentially methylated bases...")
myDiff10p = getMethylDiff(myDiff, difference=10, qvalue=0.01)

# Save the results
print("Saving results...")
write.csv(getData(myDiff), file=file.path(output_dir, "methylKit_all_diff_meth_customParams.csv"), row.names=FALSE)
write.csv(getData(myDiff10p), file=file.path(output_dir, "methylKit_sig_diff_meth_10p_customParams.csv"), row.names=FALSE)

saveRDS(list(myobj=myobj, meth=meth, myDiff=myDiff, myDiff10p=myDiff10p), 
        file=file.path(output_dir, "methylKit_unite.rds"))

print("Analysis complete.")
