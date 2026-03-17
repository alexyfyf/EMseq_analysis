# CpG Context Coverage Summary

Coverage was summarized across four CpG annotation classes from `annotatr` for `hg38`: CpG islands, shores, shelves, and open sea. Per-locus coverage was calculated from the Bismark coverage file as methylated reads plus unmethylated reads, and loci were assigned to CpG context by overlap with the annotation sets.

At the individual CpG level, the coverage boxplot showed that CpG islands did not have the highest per-site read depth relative to the other CpG contexts. In contrast, the `genomation` meta-profile, generated from regions centered and extended by +/- 1 kb, showed a clear enrichment of coverage around CpG island regions. These two observations are consistent because they describe different properties of the data: the boxplot reflects the distribution of read depth at single CpG loci, whereas the meta-profile reflects the spatial concentration of covered loci across annotated regions. Together, the results indicate that CpG islands are enriched as structured regions, even though individual CpGs within islands are not necessarily the highest-covered sites.

## Output Files

- Per-CpG coverage: [CB_G_1_per_cpg_coverage.tsv.gz](/Users/yan.a/github/emseq_coverage/CB_G_1_per_cpg_coverage.tsv.gz)
- Per-region mean coverage: [CB_G_1_per_region_mean_coverage.tsv](/Users/yan.a/github/emseq_coverage/CB_G_1_per_region_mean_coverage.tsv)
- CpG category summary: [CB_G_1_cpg_region_coverage_summary.tsv](/Users/yan.a/github/emseq_coverage/CB_G_1_cpg_region_coverage_summary.tsv)
- Plotting counts: [CB_G_1_cpg_region_plotting_counts.tsv](/Users/yan.a/github/emseq_coverage/CB_G_1_cpg_region_plotting_counts.tsv)
