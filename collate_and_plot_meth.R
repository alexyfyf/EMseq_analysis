# collate_and_plot_meth.R

library(dplyr)
library(ggplot2)
library(UpSetR)
library(GGally)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
    output_dir <- args[1]
} else {
    output_dir <- "/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq/data/test"
}

if (!dir.exists(output_dir)) {
    stop(paste("Output directory does not exist:", output_dir))
}

# 1. Load Data
message("Loading DML results...")

load_dml <- function(file_path, method_name, id_col, diff_col, is_fraction=TRUE) {
  if (file.exists(file_path)) {
    if (grepl("\\.bed$", file_path)) {
      data <- read.table(file_path, header=FALSE, sep="\t", stringsAsFactors=FALSE)
      colnames(data) <- c("chr", "start", "stop", "qvalue", "meth_diff", "n_cpgs", "pvalue", "p2d", "mean_g1", "mean_g2")
      data$id <- paste(data$chr, data$stop, sep="_")
    } else {
      data <- read.csv(file_path)
      data$id <- paste(data[, "chr"], data[, id_col], sep="_")
    }
    data$method <- method_name
    
    # Standardize meth_diff to percentage
    if (method_name == "MethylKit") {
        data$meth_diff_standardized <- data$meth.diff
    } else if (method_name == "Metilene") {
        data$meth_diff_standardized <- data$meth_diff * 100
    } else {
        data$meth_diff_standardized <- data[, diff_col] * 100
    }
    data$meth_diff_abs <- abs(data$meth_diff_standardized)
    return(data)
  }
  return(NULL)
}

mk_data <- load_dml(file.path(output_dir, "methylKit_sig_diff_meth_10p_customParams.csv"), "MethylKit", "start", "meth.diff", FALSE)
bs_data <- load_dml(file.path(output_dir, "bsmooth_sig_diff_meth_10p_loci.csv"), "BSmooth", "pos", "meth_diff")
dss_data <- load_dml(file.path(output_dir, "dss_sig_diff_meth_10p_loci.csv"), "DSS", "pos", "diff")
met_data <- load_dml(file.path(output_dir, "metilene_dmc_results_10p.bed"), "Metilene", "stop", "meth_diff")

# Check what we have
available_methods <- list()
if (!is.null(mk_data)) available_methods$MethylKit <- mk_data
if (!is.null(bs_data)) available_methods$BSmooth <- bs_data
if (!is.null(dss_data)) available_methods$DSS <- dss_data
if (!is.null(met_data)) available_methods$Metilene <- met_data

if (length(available_methods) == 0) {
  stop("No valid DML result files found in the output directory.")
}

message(paste("Successfully loaded data for:", paste(names(available_methods), collapse=", ")))

# 2. Stacked Bar Plot of Methylation Difference Categories
message("Generating stacked bar plot...")
# Function to categorize
categorize_diff <- function(diff_val) {
  ifelse(diff_val < 20, "10-20%",
         ifelse(diff_val < 50, "20-50%", "50-100%"))
}

# Create a combined dataframe for plotting the distribution
all_diffs <- do.call(bind_rows, lapply(names(available_methods), function(m) {
  data.frame(Method = m, Diff = available_methods[[m]]$meth_diff_abs)
}))
all_diffs$Category <- factor(categorize_diff(all_diffs$Diff), levels = c("10-20%", "20-50%", "50-100%"))

# Pre-calculate counts and percentages for labels
summary_df <- all_diffs %>%
  group_by(Method, Category) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Method) %>%
  mutate(Total = sum(Count),
         Percentage = Count / Total)

# Plot 1: Raw Numbers
bar_count <- ggplot(summary_df, aes(x = Method, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 3) +
  theme_minimal() +
  labs(title = "Number of DMLs by Methylation Difference", x = "Method", y = "Count of DMLs") +
  scale_fill_brewer(palette = "Set1")

# Plot 2: Percentages
bar_pct <- ggplot(summary_df, aes(x = Method, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = scales::percent(Percentage, accuracy=0.1)), position = position_fill(vjust = 0.5), size = 3) +
  theme_minimal() +
  labs(title = "Proportion of DMLs by Methylation Difference", x = "Method", y = "Percentage of DMLs") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set1")

pdf(file.path(output_dir, "dml_difference_barplot.pdf"), width = 10, height = 5)
print(bar_count)
print(bar_pct)
dev.off()

# 3. Upset Plot for Overlaps
message("Generating UpSet plot...")
listInput <- lapply(available_methods, function(df) df$id)

pdf(file.path(output_dir, "dml_overlap_upset.pdf"), width = 7, height = 5)
upset(fromList(listInput), order.by = "freq", main.bar.color = "black", sets.bar.color = "blue")
dev.off()

# 4. Pairwise Scatterplots (Correlation)
message("Generating pairwise scatterplots...")
# First, create minimalist dataframes mapping ID to the values
available_subs <- list()
if (!is.null(mk_data)) {
  sub <- mk_data[, c("id", "meth.diff")]
  colnames(sub)[2] <- "MethylKit"
  available_subs[[length(available_subs)+1]] <- sub
}
if (!is.null(bs_data)) {
  sub <- bs_data[, c("id", "meth_diff")]
  sub$BSmooth <- sub$meth_diff * 100
  available_subs[[length(available_subs)+1]] <- sub[, c("id", "BSmooth")]
}
if (!is.null(dss_data)) {
  sub <- dss_data[, c("id", "diff")]
  sub$DSS <- sub$diff * 100
  available_subs[[length(available_subs)+1]] <- sub[, c("id", "DSS")]
}
if (!is.null(met_data)) {
  sub <- met_data[, c("id", "meth_diff")]
  sub$Metilene <- sub$meth_diff * 100
  available_subs[[length(available_subs)+1]] <- sub[, c("id", "Metilene")]
}

if (length(available_subs) >= 2) {
  merged_diffs <- available_subs[[1]]
  for (i in 2:length(available_subs)) {
    merged_diffs <- inner_join(merged_diffs, available_subs[[i]], by = "id")
  }
} else {
  merged_diffs <- data.frame()
}

if(nrow(merged_diffs) > 10) {
  # Sample down if it's too large so plotting doesn't take forever or crash memory
  plot_df <- merged_diffs
  if(nrow(plot_df) > 10000) {
    message("Subsampling 10,000 common loci for faster scatterplot generation...")
    plot_df <- plot_df[sample(nrow(plot_df), 10000), ]
  }
  
  pdf(file.path(output_dir, "dml_correlation_scatterplot.pdf"), width=10, height=10)
  p <- ggpairs(plot_df[, -1], title="Correlation of Methylation Difference (%)",
               lower = list(continuous = wrap("points", alpha = 0.3, size=0.5)))
  print(p)
  dev.off()
} else {
  message("Warning: Not enough shared CpGs (< 10) or not enough methods (>= 2) found to draw a meaningful scatterplot.")
}

# 5. Extracting P-values/Q-values for comparison (MethylKit vs DSS vs Metilene)
message("Generating p-value / q-value comparisons...")

plot_pvals_dynamic <- function(available_dfs, col_mappings, filename, title) {
  # available_dfs: list of dataframes by method name
  # col_mappings: named vector mapping method name to its p/q-value column
  
  merged <- NULL
  valid_methods <- intersect(names(available_dfs), names(col_mappings))
  
  if (length(valid_methods) < 2) {
    message(paste("Warning: Not enough methods (need >= 2) providing significance values to draw", title))
    return()
  }
  
  for (m in valid_methods) {
    df_sub <- available_dfs[[m]][, c("id", col_mappings[m])]
    colnames(df_sub)[2] <- m
    if (is.null(merged)) {
      merged <- df_sub
    } else {
      merged <- inner_join(merged, df_sub, by = "id")
    }
  }
  
  if(!is.null(merged) && nrow(merged) > 10) {
    plot_data <- merged
    # Convert to -log10
    for (m in valid_methods) {
      plot_data[[m]] <- -log10(plot_data[[m]])
    }
    
    # Replace Inf with a max cap for plotting if val was 0
    vals_only <- as.matrix(plot_data[, valid_methods])
    max_log <- max(vals_only[is.finite(vals_only)], na.rm=TRUE)
    for (m in valid_methods) {
      plot_data[[m]][!is.finite(plot_data[[m]])] <- max_log + 10
    }

    if(nrow(plot_data) > 10000) {
      plot_data <- plot_data[sample(nrow(plot_data), 10000), ]
    }

    pdf(file.path(output_dir, filename), width=8, height=8)
    p2 <- ggpairs(plot_data[, valid_methods], title=title,
                  lower = list(continuous = wrap("points", alpha = 0.3, size=0.5)))
    print(p2)
    dev.off()
  } else {
    message(paste("Warning: Not enough shared CpGs (< 10) across significance-providing methods to draw", title))
  }
}

qval_cols <- c(MethylKit="qvalue", DSS="fdr", Metilene="qvalue")
pval_cols <- c(MethylKit="pvalue", DSS="pval", Metilene="pvalue")

plot_pvals_dynamic(available_methods, qval_cols, "dml_qvalue_correlation.pdf", "Correlation of Significance (-log10 Q-Value/FDR)")
plot_pvals_dynamic(available_methods, pval_cols, "dml_pvalue_correlation.pdf", "Correlation of Significance (-log10 P-Value)")

message("All available plots generated successfully in: ", output_dir)
