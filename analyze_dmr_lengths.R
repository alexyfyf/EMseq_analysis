library(tidyverse)
library(ggplot2)

# Set base directory
base_dir <- "/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq"
output_dir <- file.path(base_dir, "output")
results_file <- file.path(output_dir, "dmr_length_analysis.pdf")

# Define cell lines and comparison patterns
cell_lines <- c("CB", "KASU", "MDSL", "MOLM", "Pt")
comparisons <- c("A_vs_G", "A_vs_V", "G_vs_V")

all_data <- list()

for (cl in cell_lines) {
  for (comp in comparisons) {
    # Data directory (handling folder naming convention)
    folder_name <- paste0(cl, "_", comp)
    file_path <- file.path(output_dir, folder_name, "dss_sig_diff_meth_10p_regions.csv")
    
    if (file.exists(file_path)) {
      df <- read.csv(file_path)
      if (nrow(df) > 0) {
        df$cell_line <- cl
        df$comparison <- comp
        df$length <- df$end - df$start + 1
        df$status <- ifelse(df$diff > 0, "Hyper", "Hypo")
        all_data[[length(all_data) + 1]] <- df[, c("cell_line", "comparison", "length", "status")]
      }
    } else {
      message(paste("Warning: File not found:", file_path))
    }
  }
}

if (length(all_data) == 0) {
  stop("No DMR data found in any of the expected directories.")
}

combined_df <- bind_rows(all_data)

# Factor levels for consistent plotting
combined_df$cell_line <- factor(combined_df$cell_line, levels = cell_lines)
combined_df$status <- factor(combined_df$status, levels = c("Hyper", "Hypo"))

# Plotting
pdf(results_file, width = 15, height = 10)

# 1. Boxplot of lengths
p1 <- ggplot(combined_df, aes(x = comparison, y = length, fill = status)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_y_log10() +
  facet_grid(status ~ cell_line, scales = "free_x") +
  theme_bw() +
  labs(title = "DMR Length Distribution across Comparisons",
       subtitle = "Log10 scale for length. Rows: Methylation Status, Columns: Cell Line",
       x = "Comparison",
       y = "Length (bp)",
       fill = "Methylation Status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(angle = 0))

print(p1)

# 2. Density plot of lengths
p2 <- ggplot(combined_df, aes(x = length, color = comparison)) +
  geom_density() +
  scale_x_log10() +
  facet_grid(status ~ cell_line) +
  theme_bw() +
  labs(title = "DMR Length Density",
       subtitle = "Rows: Methylation Status, Columns: Cell Line",
       x = "Length (bp) (Log10 scale)",
       y = "Density") +
  theme(strip.text.y = element_text(angle = 0))

print(p2)

dev.off()

# Save summary stats
summary_stats <- combined_df %>%
  group_by(cell_line, comparison, status) %>%
  summarise(
    count = n(),
    mean_length = mean(length),
    median_length = median(length),
    sd_length = sd(length),
    min_length = min(length),
    max_length = max(length),
    .groups = "drop"
  )

write.csv(summary_stats, file.path(output_dir, "dmr_length_summary.csv"), row.names = FALSE)

message(paste("Analysis complete. Plots saved to:", results_file))
message(paste("Summary stats saved to:", file.path(output_dir, "dmr_length_summary.csv")))
