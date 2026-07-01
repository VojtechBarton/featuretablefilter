# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Multi-type input support for `run_filtering_pipeline()`**: Pipeline now accepts `data.frame`, `matrix`, `phyloseq`, or `TreeSummarizedExperiment` objects directly (not just file paths)
- **Automatic class preservation**: Output `filtered_table` returns in the same class as input, preserving metadata (taxonomy, phylogenetic tree, sample data for phyloseq; rowData, colData, reducedDims for TSE)
- **New `input_class` return value**: Pipeline result includes the class of the input object for programmatic access
- **`SECURITY.md`**: Security vulnerability reporting policy
- **`CITATION.cff`**: Standard citation file for software citation
- **GitHub Actions workflows**:
  - `R-CMD-check.yaml`: Cross-platform R CMD check (macOS, Windows, Ubuntu)
  - `test-coverage.yaml`: Code coverage reporting via Codecov
  - `pkgdown.yaml`: Automatic documentation site deployment
  - `docker.yaml`: Docker image build and push to GHCR

### Changed

- Renamed `input_file` parameter to `input` in `run_filtering_pipeline()` to support both file paths and objects
- Updated all documentation examples from `feature-table.tsv` to `example_feature_table.tsv`
- Removed sensitive production data (`feature-table.tsv`) from repository

### Fixed

- Bug in report generation where `input_file` variable was referenced after renaming to `input`

---

## [0.1.0] - 2026-06-15

### Added

#### Core Package Infrastructure
- Initial package structure with `DESCRIPTION`, `NAMESPACE`, `LICENSE` (GPL-3)
- Complete README with installation instructions, quick start examples, and feature overview
- CHANGELOG.md for tracking changes
- CLAUDE.md for AI assistant context
- Docker configuration (`Dockerfile`, `docker-compose.yml`) for reproducible analysis environment
- Vignette: `featuretablefilter.Rmd` with comprehensive usage guide

#### Core Filtering Functions
- **`filter_by_coverage()`**: Remove low-coverage samples using absolute minimum read thresholds
- **`filter_features_by_abundance()`**: Filter features by absolute or relative abundance with `min_samples` support
- **`filter_by_relative_cutoff()`**: Relative threshold filtering based on minimum-coverage sample

#### Coverage Estimation Functions
- **`estimate_mad_cutoff()`**: MAD-based threshold estimation (median ± multiplier × MAD)
- **`estimate_iqr_cutoff()`**: IQR/Tukey-based threshold estimation (fences method)
- **`estimate_good_coverage()`**: Good's coverage estimator for ecological completeness
- **`estimate_chao_coverage()`**: Chao's coverage estimator (conservative completeness)
- **`filter_by_coverage_estimator()`**: Filter samples by ecological coverage threshold (Good's or Chao's)

#### QC & Metrics
- **`compute_filtering_qc()`**: Comprehensive QC metrics including:
  - Sparsity changes (original vs filtered)
  - Retention rates (reads, features, samples)
  - Rank-abundance stability (top N overlap, Jaccard similarity, Spearman correlation)
  - Procrustes analysis for compositional similarity (via vegan)
- **`calculate_sparsity()`**: Calculate sparsity metric (percentage of zeros)
- **`calculate_cv()`**: Calculate coefficient of variation per feature

#### Visualization Functions
- **`plot_qc_comparison()`**: Generate 7 comparison plots:
  1. Sample coverage distribution (log10, faceted)
  2. Feature abundance distribution (log10, faceted)
  3. Per-sample sparsity histogram
  4. Retention rates barplot
  5. Top features stacked barplot
  6. Heatmaps (original and filtered)
  7. Presence frequency comparison
- **`plot_presence_analysis()`**: Feature prevalence and sample richness histograms
- **`plot_coverage_histogram()`**: Sample coverage histograms with optional threshold line
- **`plot_top_features_stacked()`**: Top N relative abundance stacked bar charts
- **`plot_scree()`**: Scree/saturation diagnostic plots for filtering thresholds

#### Pipeline & Utilities
- **`run_filtering_pipeline()`**: Complete filtering workflow with:
  - Configurable methods per step (coverage, abundance)
  - File outputs (filtered table, plots, text report)
  - Progress verbosity control
- **`load_feature_table()`**: Auto-detects TSV/CSV format, header detection, feature column identification
- **`convert_to_relative()`**: Convert count data to relative abundances

---

## [0.1.1] - 2026-06-18

### Added

- **`calc_shannon_ens()`**: Shannon diversity converted to effective species count (Hill number q=1)
- **`calc_simpson_ens()`**: Simpson diversity as effective species count (Hill number q=2)
- **`calc_hill_numbers()`**: Full diversity profile across multiple q values
- Enhanced `compute_filtering_qc()` to include ENS (Effective Number of Species) retention percentages

### Changed

- Updated `compute_filtering_qc()` to report diversity retention alongside traditional metrics

---

## [0.1.2] - 2026-06-22

### Added

- **`identify_sparsity_elbow()`**: Detect elbow point in richness-depth curve using:
  - Kneedle algorithm
  - Maximum derivative method
  - Second derivative method
- **`plot_sparsity_elbow()`**: Visualize richness-depth relationship with elbow marker
- **`analyze_depth_sparsity()`**: Full depth-sparsity analysis with MAD/IQR outlier detection
- **`plot_depth_sparsity()`**: Visualize depth-sparsity relationship with outliers flagged
- **`filter_depth_sparsity_outliers()`**: Convenience function to remove outlier samples
- **`plot_reads_vs_asvs()`**: Total Reads vs Observed ASVs plot with MAD outlier flagging

### Changed

- Improved `plot_presence_analysis()` with proper histogram bins instead of barplots
- Fixed feature prevalence / sample richness calculation (swapped rowSums/colSums)

### Fixed

- Y-axis count display issues in presence analysis plots
- Missing legend entries for mean/median statistics
- Duplicate vertical lines in facet panels

---

## [0.1.3] - 2026-06-25

### Added

- **`compute_mutual_information()`**: Pairwise mutual information using k-nearest neighbor (KSG) estimator
- **`analyze_feature_network()`**: Network centrality metrics (degree, strength, betweenness)
- **`filter_by_network_connectivity()`**: Remove features with no significant connections (likely artifacts)
- **`plot_feature_network()`**: Visualize feature connectivity network using igraph
- **`filter_cross_talk()`**: Index hopping / cross-talk detection and filtering with modes:
  - `zero`: Set suspected leakage reads to zero
  - `remove_feature`: Remove entire feature if any leakage detected
  - `flag`: Flag leakage but don't modify data

### Added

- **`filter_by_singleton_ratio()`**: PCR artifact detection by filtering samples with high singleton/doubleton ratios

---

## [0.1.4] - 2026-06-26

### Added

- **`filter_features_joint()`**: Joint abundance-prevalence filtering with AND/OR logic:
  - OR: Keep features meeting EITHER abundance OR prevalence criteria
  - AND: Keep features meeting BOTH abundance AND prevalence criteria
  - Detailed breakdown of retention reasons (abundance only, prevalence only, both, neither)

### Changed

- Enhanced `run_filtering_pipeline()` to support joint filtering method

---

## [0.1.5] - 2026-06-29

### Added

- **`from_phyloseq()`**: Convert phyloseq object to data.frame (with optional taxonomy inclusion)
- **`to_phyloseq()`**: Create phyloseq from data.frame (with optional taxonomy, tree, sample data)
- **`from_TSE()`**: Convert TreeSummarizedExperiment/SingleCellExperiment to data.frame
- **`to_TSE()`**: Create TreeSummarizedExperiment from data.frame (with rowData, colData, reducedDims, rowTree, rowLinks)
- **`convert_feature_table()`**: Generic auto-detect conversion between any supported format

### Changed

- All functions now gracefully handle optional dependencies (phyloseq, TreeSummarizedExperiment, vegan) with warnings

---

## [0.1.6] - 2026-07-01

### Added

- **`example_feature_table.tsv`**: Synthetic example dataset (50 features × 20 samples) for documentation and testing

### Changed

- Removed sensitive production data (`feature-table.tsv`) from repository
- Updated all documentation examples to use `example_feature_table.tsv`

### Fixed

- Bug in report generation where `input_file` variable was referenced after renaming

---

## Development History Summary

| Version | Date | Key Additions |
|---------|------|---------------|
| 0.1.0 | 2026-06-15 | Initial release: core filtering, coverage estimation, QC metrics, visualization, pipeline |
| 0.1.1 | 2026-06-18 | Hill numbers / Effective Species Count diversity metrics |
| 0.1.2 | 2026-06-22 | Sparsity elbow detection, depth-sparsity outlier analysis |
| 0.1.3 | 2026-06-25 | Network-based filtering, mutual information, cross-talk filtering, singleton ratio filtering |
| 0.1.4 | 2026-06-26 | Joint abundance-prevalence filtering with AND/OR logic |
| 0.1.5 | 2026-06-29 | Native phyloseq and TreeSummarizedExperiment support |
| 0.1.6 | 2026-07-01 | Example data, GitHub standards, multi-type input support |

---

[Unreleased]: https://github.com/VojtechBarton/featuretablefilter/compare/v0.1.6...HEAD
[0.1.0]: https://github.com/VojtechBarton/featuretablefilter/releases/tag/v0.1.0
