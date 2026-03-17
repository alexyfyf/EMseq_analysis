# annotate_results.R
# This script acts as a wrapper to call the user-provided annotation function.

library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(annotatr)

# --- USER PROVIDED CORE FUNCTION ---
# Sourced from: /home/users/allstaff/yan.a/yan.a/software/utils/rrbs/comp_meth_anno.R
template_path <- "/home/users/allstaff/yan.a/yan.a/software/utils/rrbs/comp_meth_anno.R"
if (file.exists(template_path)) {
    message(paste("Sourcing core annotation function from:", template_path))
    source(template_path)
} else {
    stop(paste("Template script not found at:", template_path))
}

# --- SCRIPT LOGIC ---

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
    output_dir <- args[1]
} else {
    stop("Usage: Rscript annotate_results.R <output_dir>")
}

if (!dir.exists(output_dir)) {
    stop(paste("Output directory does not exist:", output_dir))
}

genome <- "hg38"
message(paste("Starting annotation for directory:", output_dir))

# Build CpG annotations once for performance
message("Building CpG annotations...")
cpgs_info <- build_annotations(genome = genome, annotations = paste0(genome, "_cpgs"))

# Helper function to convert data frame to GRanges
df_to_gr <- function(df, seqnames_col="chr", start_col="start", end_col="end") {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  makeGRangesFromDataFrame(df, 
                           keep.extra.columns = TRUE, 
                           seqnames.field = seqnames_col, 
                           start.field = start_col, 
                           end.field = end_col)
}

# Define result files to process with their specific column mappings
dml_files <- list(
  list("methylKit_sig_diff_meth_10p_customParams.csv", "MethylKit", "chr", "start", "end", "meth.diff"),
  list("bsmooth_sig_diff_meth_10p_loci.csv", "BSmooth", "chr", "pos", "pos", "meth_diff"),
  list("dss_sig_diff_meth_10p_loci.csv", "DSS", "chr", "pos", "pos", "diff"),
  list("metilene_dmc_results_10p.bed", "Metilene", "V1", "V2", "V3", "V5")
)

dmr_files <- list(
  list("dss_sig_diff_meth_10p_regions.csv", "DSS_DMR", "chr", "start", "end", "diff.Methy"),
  list("metilene_dmr_results_10p.bed", "Metilene_DMR", "V1", "V2", "V3", "V5")
)

process_and_annotate <- function(file_info, is_dml=TRUE) {
  file_path <- file.path(output_dir, file_info[[1]])
  label <- file_info[[2]]
  
  if (!file.exists(file_path)) {
    message(paste("Skipping", label, "- file not found:", file_info[[1]]))
    return(NULL)
  }
  
  message(paste("Processing", label, "..."))
  
  if (grepl("\\.csv$", file_path)) {
    data <- read.csv(file_path)
  } else {
    data <- read.table(file_path, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  }
  
  if (nrow(data) == 0) {
    message(paste("Skipping", label, "- file is empty."))
    return(NULL)
  }
  
  diff_col <- file_info[[6]]
  
  # Check if diff column exists
  if (!(diff_col %in% colnames(data))) {
     message(paste("Warning: Column", diff_col, "not found in", label, "- available columns are:", paste(colnames(data), collapse=", ")))
     # Try to fallback or skip
     return(NULL)
  }

  gr <- df_to_gr(data, file_info[[3]], file_info[[4]], file_info[[5]])
  
  if (!is.null(gr)) {
    message(paste("Calling user function 'comp_meth_anno' for", label, "..."))
    
    # CALLING THE PROVIDED FUNCTION with tryCatch
    tryCatch({
      anno_res <- comp_meth_anno(gr, 
                                 cpgs_info = cpgs_info, 
                                 column = diff_col, 
                                 simplify = FALSE, 
                                 plot = "bar", 
                                 genome = genome)
      
      # Save outputs to the same folder
      prefix <- if(is_dml) "dml" else "dmr"
      out_prefix <- file.path(output_dir, paste0(prefix, "_annotated_", tolower(label)))
      
      pdf(paste0(out_prefix, "_plots.pdf"), width=10, height=12)
      print(anno_res$plot)
      dev.off()
      
      write.csv(anno_res$summary, paste0(out_prefix, "_summary.csv"))
      write.csv(as.data.frame(anno_res$tss), paste0(out_prefix, "_details.csv"), row.names = FALSE)
      
      message(paste("Successfully annotated and saved results for", label))
    }, error = function(e) {
      message(paste("Error annotating", label, ":", e$message))
    })
  }
}

# Execute annotations
message("\n--- Annotating DMLs ---")
for (f in dml_files) {
  process_and_annotate(f, is_dml=TRUE)
}

message("\n--- Annotating DMRs ---")
for (f in dmr_files) {
  process_and_annotate(f, is_dml=FALSE)
}

message("\nAll available annotations completed.")
