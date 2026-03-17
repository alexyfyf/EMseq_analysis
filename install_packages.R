# install_packages.R
# Script to install missing Bioconductor packages for R 4.4.2

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")

packages <- c("annotatr", "ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene")

# Identify which ones are really missing
missing <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]

if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    BiocManager::install(missing, ask = FALSE, update = FALSE)
} else {
    message("All required packages are already installed.")
}
