# featuretablefilter NEWS

## featuretablefilter 1.0.0 (Unreleased)

### Added

- **Multi-type input support for `run_filtering_pipeline()`**: Pipeline now accepts `data.frame`, `matrix`, `phyloseq`, or `TreeSummarizedExperiment` objects directly (not just file paths)
- **Automatic class preservation**: Output `filtered_table` returns in the same class as input, preserving metadata (taxonomy, phylogenetic tree, sample data for phyloseq; rowData, colData, reducedDims for TSE)
- **New `input_class` return value**: Pipeline result includes the class of the input object for programmatic access
- Network-based filtering functions:
  - `compute_mutual_information()`: Pairwise MI using k-nearest neighbor (KSG) estimator
  - `analyze_feature_network()`: Network centrality metrics (degree, strength, betweenness)
  - `filter_by_network_connectivity()`: Remove disconnected features (likely artifacts)
  - `plot_feature_network()`: Visualize feature connectivity network
- Sparsity elbow detection functions:
  - `identify_sparsity_elbow()`: Detect elbow point using kneedle, max_derivative, or second_derivative methods
  - `plot_sparsity_elbow()`: Visualize richness-depth relationship with elbow marker
  - `analyze_depth_sparsity()`: Full depth-sparsity analysis with MAD/IQR outlier detection
  - `plot_depth_sparsity()`: Visualize depth-sparsity relationship
  - `filter_depth_sparsity_outliers()`: Remove depth-sparsity outliers
  - `plot_reads_vs_asvs()`: Total Reads vs Observed ASVs plot with outlier flagging
- Diversity retention metrics:
  - `calc_shannon_ens()`: Shannon diversity as effective species count (Hill q=1)
  - `calc_simpson_ens()`: Simpson diversity as effective species count (Hill q=2)
  - `calc_hill_numbers()`: Full diversity profile across multiple q values
- Native phyloseq and TreeSummarizedExperiment support:
  - `from_phyloseq()` / `to_phyloseq()`: Convert between phyloseq and data.frame
  - `from_TSE()` / `to_TSE()`: Convert between TreeSummarizedExperiment and data.frame
  - `convert_feature_table()`: Generic auto-detect conversion between any supported format
- Cross-talk/index-hopping filtering:
  - `filter_cross_talk()`: Remove suspected index-hopping artifacts
  - `filter_index_hopping()`: Alias for filter_cross_talk
- Joint abundance-prevalence filtering:
  - `filter_features_joint()`: Combined filtering with AND/OR logic
- Coverage estimation functions:
  - `estimate_good_coverage()`: Good's coverage estimator
  - `estimate_chao_coverage()`: Chao's coverage estimator
  - `filter_by_coverage_estimator()`: Filter by ecological coverage threshold
- Scree/saturation analysis:
  - `compute_scree()`: Evaluate retention across threshold gradients
  - `plot_scree()`: Visualize scree diagnostics
- Interactive Shiny dashboard:
  - `runDashboard()`: Launch interactive web-based GUI for exploratory filtering
- Comprehensive QC metrics:
  - `compute_filtering_qc()`: Sparsity, retention rates, rank-abundance stability, Procrustes analysis, ENS diversity retention
- Visualization functions:
  - `plot_qc_comparison()`: 7 comparison plots (coverage, abundance, sparsity, retention, heatmaps, stacked bars)
  - `plot_presence_analysis()`: Feature prevalence and sample richness histograms
  - `plot_coverage_histogram()`: Sample coverage histograms with threshold line
  - `plot_top_features_stacked()`: Top N relative abundance stacked bars
  - `plot_scree()`: Scree/saturation diagnostic plots

### Changed

- Renamed `input_file` parameter to `input` in `run_filtering_pipeline()` to support both file paths and objects
- Updated all documentation examples from `feature-table.tsv` to `example_feature_table.tsv`

### Fixed

- Bug in report generation where `input_file` variable was referenced after renaming to `input`

---

## featuretablefilter 0.1.0 (2025-06-15)

### Added

- Initial package release with core filtering functionality
- Basic vignette and documentation
- Example feature table for testing and demonstration
