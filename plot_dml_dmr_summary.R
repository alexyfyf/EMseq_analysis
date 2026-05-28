library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

output_dir <- "/home/users/allstaff/yan.a/davidson_longread/yan.a/XZ/emseq/output"
comparisons <- list.dirs(output_dir, recursive = FALSE, full.names = FALSE)
# Filter out non-comparison folders
comparisons <- comparisons[!comparisons %in% c("compare_dss", "cpg_coverage_plots") & !grepl("\\.pdf|\\.csv", comparisons)]

results <- data.frame()

# Function to efficiently get binned counts
get_binned_counts <- function(file_path, col_idx, has_header = TRUE, is_pct = FALSE) {
  if (file.exists(file_path)) {
    res <- tryCatch({
      # fread is very fast for reading specific columns
      dt <- fread(file_path, select = col_idx, header = has_header, data.table = FALSE)
      if (nrow(dt) == 0) {
        return(c("0-20%" = 0, "20-50%" = 0, "50-100%" = 0))
      }
      
      # Extract column and compute absolute difference
      vals <- abs(as.numeric(dt[[1]]))
      if (is_pct) vals <- vals / 100
      
      bin1 <- sum(vals <= 0.20, na.rm = TRUE)
      bin2 <- sum(vals > 0.20 & vals <= 0.50, na.rm = TRUE)
      bin3 <- sum(vals > 0.50, na.rm = TRUE)
      
      c("0-20%" = bin1, "20-50%" = bin2, "50-100%" = bin3)
    }, error = function(e) {
      c("0-20%" = 0, "20-50%" = 0, "50-100%" = 0)
    })
    return(res)
  } else {
    return(c("0-20%" = 0, "20-50%" = 0, "50-100%" = 0))
  }
}

cat("Processing comparisons...\n")

for (comp in comparisons) {
  comp_dir <- file.path(output_dir, comp)
  
  # DML binned counts
  bsmooth_dml <- get_binned_counts(file.path(comp_dir, "bsmooth_sig_diff_meth_10p_loci.csv"), 5, TRUE, FALSE)
  dss_dml <- get_binned_counts(file.path(comp_dir, "dss_sig_diff_meth_10p_loci.csv"), 5, TRUE, FALSE)
  methylkit_dml <- get_binned_counts(file.path(comp_dir, "methylKit_sig_diff_meth_10p_customParams.csv"), 7, TRUE, TRUE)
  metilene_dml <- get_binned_counts(file.path(comp_dir, "metilene_dmc_results_10p.bed"), 5, FALSE, FALSE)
  
  # DMR binned counts
  bsmooth_dmr <- get_binned_counts(file.path(comp_dir, "bsmooth_sig_diff_meth_10p_regions.csv"), 5, TRUE, FALSE)
  dss_dmr <- get_binned_counts(file.path(comp_dir, "dss_sig_diff_meth_10p_regions.csv"), 8, TRUE, FALSE)
  metilene_dmr <- get_binned_counts(file.path(comp_dir, "metilene_dmr_results_10p.bed"), 5, FALSE, FALSE)
  methylkit_dmr <- c("0-20%" = 0, "20-50%" = 0, "50-100%" = 0) # methylKit doesn't output DMRs in this format
  
  methods <- c("BSmooth", "DSS", "methylKit", "metilene")
  
  # Bind DML
  dml_mat <- rbind(bsmooth_dml, dss_dml, methylkit_dml, metilene_dml)
  for(i in 1:length(methods)) {
    results <- rbind(results, data.frame(
      Comparison = comp,
      Method = methods[i],
      Feature = "DML",
      MethDiff = c("0-20%", "20-50%", "50-100%"),
      Count = dml_mat[i, ]
    ))
  }
  
  # Bind DMR
  dmr_mat <- rbind(bsmooth_dmr, dss_dmr, methylkit_dmr, metilene_dmr)
  for(i in 1:length(methods)) {
    results <- rbind(results, data.frame(
      Comparison = comp,
      Method = methods[i],
      Feature = "DMR",
      MethDiff = c("0-20%", "20-50%", "50-100%"),
      Count = dmr_mat[i, ]
    ))
  }
}

cat("Done processing. Plotting results...\n")

# Prepare data for plotting
results <- results %>%
  mutate(
    Method = factor(Method, levels = c("BSmooth", "DSS", "methylKit", "metilene")),
    MethDiff = factor(MethDiff, levels = c("50-100%", "20-50%", "0-20%")), # Stack order
    CellLine = sub("_.*", "", Comparison)
  )

# Create plot
p <- ggplot(results, aes(x = Method, y = Count, fill = MethDiff)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  facet_grid(Feature ~ CellLine, scales = "free") +
  scale_fill_manual(values = c("50-100%" = "#e41a1c", "20-50%" = "#377eb8", "0-20%" = "#4daf4a")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  ) +
  labs(title = "Number of DML and DMR by Method across Comparisons", 
       subtitle = "Broken down by Absolute Methylation Difference",
       x = "Method", 
       y = "Count",
       fill = "Meth Diff")

# Save outputs
data_dir <- "/home/users/allstaff/yan.a/davidson_longread/yan.a/XZ/emseq/data"
dir.create(data_dir, showWarnings = FALSE)

plot_file <- file.path(data_dir, "dml_dmr_summary_plot.pdf")
csv_file <- file.path(data_dir, "dml_dmr_summary_counts.csv")

ggsave(plot_file, p, width = 14, height = 10)
write.csv(results, csv_file, row.names = FALSE)

cat("Successfully saved plot to:", plot_file, "\n")
cat("Successfully saved data to:", csv_file, "\n")
