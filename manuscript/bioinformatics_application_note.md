# featuretablefilter: an R/Bioconductor package for reproducible filtering of microbiome feature tables

**Vojtech Barton**

Affiliation: *Department of [Field], [Institution], [City, Country]*

Corresponding author: Vojtech Barton, vojtech.barton@gmail.com

---

## Abstract

**Motivation:** Microbiome feature tables produced by amplicon sequencing pipelines are inherently sparse, noisy and affected by technical artifacts such as low-coverage samples, rare features, PCR errors and index hopping. Although filtering is one of the most consequential preprocessing steps in microbiome analysis, existing tools are fragmented across packages and rarely provide integrated quality-control (QC) metrics or transparent reporting of filtering decisions.

**Results:** We present `featuretablefilter`, an R package distributed through Bioconductor that unifies sample coverage, feature abundance, singleton-ratio, cross-talk, depth-sparsity and network-based filtering in a single reproducible workflow. The package offers native support for `phyloseq` and `TreeSummarizedExperiment` objects, data-driven threshold estimation, comprehensive QC metrics and publication-ready visualizations. A single command, `run_filtering_pipeline()`, executes the full workflow while preserving the input object class and generating traceable reports. An interactive Shiny dashboard (`runDashboard()`) allows users to explore filtering parameters and inspect QC metrics in real time without writing code.

**Availability:** `featuretablefilter` is available from Bioconductor (https://bioconductor.org/packages/featuretablefilter) and GitHub (https://github.com/VojtechBarton/featuretablefilter) under the Artistic License 2.0. Docker images and interactive Shiny dashboard are included.

---

## 1 Introduction

Microbiome studies rely on feature tables that summarize the abundance of amplicon sequence variants (ASVs), operational taxonomic units (OTUs) or other taxonomic features across samples. These tables are typically sparse, containing many rare or zero-abundance entries, and are sensitive to sequencing depth, PCR/sequencing errors and barcode leakage (Callahan *et al.*, 2016; Stewart, 2019). Preprocessing decisions—especially which samples and features to retain—can strongly influence downstream diversity, differential abundance and network analyses (Weiss *et al.*, 2017; Nearing *et al.*, 2022).

Despite the importance of filtering, many microbiome analyses still apply thresholds in an ad hoc manner using custom scripts or isolated functions from multiple packages. This fragmentation makes it difficult to compare filtering choices, reproduce workflows, or evaluate their impact on data composition. `featuretablefilter` addresses this gap by providing a cohesive, extensible and well-documented R/Bioconductor package for filtering feature tables with integrated QC, visualization and reporting.

## 2 Implementation

`featuretablefilter` is implemented in R (≥4.1.0) and integrates with the Bioconductor ecosystem. The package imports core infrastructure including `ggplot2` for visualization and `S4Vectors` for tabular data handling, while `phyloseq`, `TreeSummarizedExperiment`, `vegan`, `ComplexHeatmap` and other packages are listed as suggested dependencies for extended functionality.

The package is organized around three layers:

1. **Atomic filtering functions** implement individual filtering strategies:
   * `filter_by_coverage()` removes low-coverage samples using absolute, MAD or IQR thresholds.
   * `filter_by_coverage_estimator()` uses ecological completeness estimators (Good's or Chao's coverage) to retain samples that are sufficiently sequenced.
   * `filter_features_by_abundance()` and `filter_features_joint()` perform absolute, relative or joint abundance–prevalence filtering with AND/OR logic.
   * `filter_by_singleton_ratio()` flags samples with excessive singletons/doubletons, a signature of PCR or sequencing artifacts.
   * `filter_cross_talk()` detects and mitigates index-hopping leakage.
   * `filter_depth_sparsity_outliers()` and `filter_by_network_connectivity()` identify outlier samples or spurious features based on richness–depth relationships and mutual-information networks, respectively.

2. **Diagnostic and QC functions** support data-driven threshold selection. `estimate_mad_cutoff()`, `estimate_iqr_cutoff()` and `identify_sparsity_elbow()` estimate cutoffs from the data itself; `compute_filtering_qc()` computes sparsity, retention rates, rank-abundance stability (Jaccard, Spearman), Hill numbers and Procrustes similarity between original and filtered tables.

3. **The integrated pipeline** `run_filtering_pipeline()` orchestrates loading, filtering, QC, visualization and reporting in one call. Inputs can be `data.frame`, TSV/CSV files, `phyloseq` or `TreeSummarizedExperiment` objects; outputs preserve the original class, including metadata and phylogenetic trees. The pipeline optionally produces a plain-text, Markdown or PDF report, coverage and abundance plots, presence–frequency analyses, scree/sensitivity analyses and a Shiny dashboard (`runDashboard()`).

A browser-based **Shiny dashboard** can be launched with `runDashboard()`. It exposes the main filtering parameters through an interactive interface and updates coverage histograms, abundance distributions, presence–frequency plots and QC comparisons as users adjust thresholds. This is particularly useful for collaborators unfamiliar with R, for teaching and for rapid exploration of the sensitivity of downstream results to filtering choices.

All functions are documented with `roxygen2`, unit-tested with `testthat` and illustrated in two package vignettes. The source repository provides a Dockerfile for reproducible deployment.

## 3 Results and Discussion

`featuretablefilter` enables systematic exploration of filtering choices rather than reliance on fixed thresholds. For example, coverage filtering can be driven by absolute read counts, robust outlier statistics (MAD/IQR) or ecological coverage estimators that account for unseen diversity. The joint abundance–prevalence filter lets users retain features that are either consistently present across samples or highly abundant when present, preserving ecological core members while reducing rare-noise features. Cross-talk and singleton-ratio filters help address well-documented technical artifacts in Illumina sequencing.

QC metrics and visualizations are generated before and after filtering, allowing researchers to assess how many reads, features and samples are retained, whether top-ranked features remain stable, and how compositional structure changes. Scree analysis systematically evaluates the effect of different thresholds, helping users avoid over-filtering. The interactive Shiny dashboard complements these programmatic tools by letting users adjust cutoffs and immediately see how retention rates, sparsity and top-feature composition respond. The package ships with synthetic example datasets (`example_feature_table`, `example_phyloseq`, `example_treesummarizedexperiment`) so users can evaluate functionality without external data.

By returning the same object class as the input, `featuretablefilter` integrates cleanly with downstream Bioconductor workflows such as differential abundance (`DESeq2`, `edgeR`), ordination (`vegan`) and network inference. The generated reports make filtering decisions transparent and auditable, supporting reproducibility requirements in microbiome research.

## 4 Conclusion

`featuretablefilter` provides a unified, reproducible and extensible framework for filtering microbiome feature tables. Its combination of multiple filtering strategies, native support for standard Bioconductor classes, integrated QC metrics and interactive reporting makes it a practical tool for both exploratory and production microbiome analyses.

---

## Availability

* Package source and documentation: https://bioconductor.org/packages/featuretablefilter
* Development repository and issue tracker: https://github.com/VojtechBarton/featuretablefilter
* Docker image: provided in the repository (`Dockerfile` and `docker-compose.yml`)
* Interactive dashboard: launch with `runDashboard()` after installing the package
* License: Artistic License 2.0

---

## Funding

*Funding information to be inserted by the author(s).*

---

## References

1. Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J., & Holmes, S. P. (2016). DADA2: high-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13(7), 581–583.

2. McMurdie, P. J., & Holmes, S. (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. *PLoS ONE*, 8(4), e61217.

3. Morgan, M., Obenchain, V., Hester, J., & Pagès, H. (2026). *Bioconductor: Open source software for genomics*. https://bioconductor.org

4. Nearing, J. T., Douglas, G. M., Comeau, A. M., & Langille, M. G. I. (2022). Denooting the effects of sequencing depth and conditioning on the analysis of microbiome data. *PeerJ*, 10, e13061.

5. Stewart, C. J. (2019). Microbiome studies need to adopt a more robust experimental design. *Microbiome*, 7, 105.

6. Weiss, S., Xu, Z. Z., Peddada, S., Amir, A., Bittinger, K., Gonzalez, A., Lozupone, C., Zaneveld, J. R., Vázquez-Baeza, Y., Birmingham, A., Hyde, E. R., & Knight, R. (2017). Normalization and microbial differential abundance strategies depend upon data characteristics. *Microbiome*, 5, 27.

