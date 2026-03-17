# EM-seq Differential Methylation Analysis Scripts

This directory contains the scripts for performing differential methylation analysis (DML/DMR) on EM-seq data. The pipeline integrates multiple tools (`MethylKit`, `BSmooth`, `DSS`, and `Metilene`) and provides annotation and visualization.

## Workflow Execution Order

To run the complete analysis, follow this order:

### 1. Pre-analysis / Quality Control (Optional)
*   **`run_plot_cpg_coverage.sh`**: (Slurm Array) Generates CpG coverage plots to verify data quality before deep analysis.

### 2. Main Differential Methylation Analysis
*   **`run_comparisons.sh`**: **(Primary Master Script)** A Slurm array job that performs pairwise treatment comparisons (e.g., A vs V, A vs G, G vs V) across all specified cell lines. 
    *   This script coordinates the execution of `MethylKit`, `BSmooth`, `DSS`, and `Metilene`.
    *   It also runs an initial visualization (`collate_and_plot_meth.R`) for each comparison.

### 3. Annotation & Functional Analysis
*   **`run_annotation_all.sh`**: Annotates the findings from the main analysis with genomic features (promoters, exons, etc.).
*   **`submit_compare_cell_line_dss.sh`**: Performs cross-comparison analysis within each cell line (e.g., comparing overlap between A vs G and A vs V) using UpSet and Sankey plots.

### 4. Summary & Statistics
*   **`analyze_dmr_lengths.R`**: Generates global summary statistics and plots for DMR lengths across the entire project.

---

## Script Descriptions

### Shell Scripts (Slurm Wrappers)
- `run_comparisons.sh`: Master array job for parallel differential analysis.
- `run_annotation_all.sh`: Batch script to annotate all comparison results.
- `run_plotting_all.sh`: Utility to regenerate plots for all comparisons without re-running analysis.
- `run_plot_cpg_coverage.sh`: Array job for CpG coverage visualization.
- `submit_compare_cell_line_dss.sh`: Wrapper for multi-treatment comparison analysis.
- `run_diff_metilene.sh`: Wrapper specifically for the `metilene` tool.

### R Scripts (Analysis & Logic)
- `run_diff_methylkit.R`: Performs DML/DMR analysis using the `MethylKit` package.
- `run_diff_bsmooth.R`: Performs analysis using `BSmooth`.
- `run_diff_dss.R`: Performs analysis using the `DSS` package.
- `collate_and_plot_meth.R`: Visualizes results from multiple tools for a single comparison.
- `annotate_results.R`: Handles genomic annotation using `ChIPpeakAnno`.
- `compare_cell_line_dss.R`: Logic for cross-comparison (3-way/2-way) overlaps and Sankey plots.
- `analyze_dmr_lengths.R`: Aggregates and plots DMR length distributions.
- `plot_cpg_region_coverage.R`: Generates coverage plots for QC.
- `install_packages.R`: Helper script to ensure all R dependencies are met.

## Logs and Output
- **Logs**: Standard output and error logs are stored in the `logs/` subdirectory.
- **Results**: Analysis results (CSV files and plots) are generally stored in the project's `output/` directory, organized by comparison name.
