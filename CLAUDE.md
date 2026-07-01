# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

**featuretablefilter** - An R package for filtering microbiome feature tables based on coverage, abundance, and network connectivity criteria. Supports `data.frame`, `phyloseq`, and `TreeSummarizedExperiment` objects natively.

## Development Commands

```bash
# Run all tests
Rscript -e "testthat::test_dir('tests/testthat')"

# Run a specific test file
Rscript -e "testthat::test_file('tests/testthat/test-compute_qc_metrics.R')"

# Run a specific test by name
Rscript -e "testthat::test_file('tests/testthat/test-compute_qc_metrics.R', filter = 'sparsity')"

# Build and install package locally
R CMD INSTALL .

# Check package
R CMD check .

# Document functions with roxygen2
Rscript -e "roxygen2::roxygenise()"

# Build package tarball
R CMD build .
```

## Architecture

The package follows a modular structure with distinct functional categories:

### Core Filtering Functions
- **`filter_coverage.R`** (`filter_by_coverage`) - Removes low-coverage samples using absolute thresholds
- **`filter_features.R`** (`filter_features_by_abundance`) - Filters features by absolute or relative abundance with `min_samples` support
- **`filter_relative_cutoff.R`** (`filter_by_relative_cutoff`) - Relative cutoff filtering based on min-coverage sample

### Coverage Estimation Functions
- **`estimate_mad_cutoff.R`** (`estimate_mad_cutoff`) - MAD-based threshold estimation (median ± multiplier × MAD)
- **`estimate_iqr_cutoff.R`** (`estimate_iqr_cutoff`) - IQR/Tukey-based threshold estimation
- **`estimate_coverage_metrics.R`** (`estimate_good_coverage`, `estimate_chao_coverage`) - Ecological completeness estimators
- **`filter_by_coverage_estimator.R`** (`filter_by_coverage_estimator`) - Filter samples by Good's or Chao's coverage threshold

### Singleton Ratio & Cross-Talk Filtering
- **`filter_by_singleton_ratio.R`** (`filter_by_singleton_ratio`) - PCR artifact detection via singleton/doubleton ratios
- **`filter_cross_talk.R`** (`filter_cross_talk`, `filter_index_hopping`) - Index hopping/cross-contamination filtering

### Network-Based Filtering
- **`filter_by_network.R`** - Information theory and graph theory based filtering
  - `compute_mutual_information()` - Pairwise MI using k-nearest neighbor (KSG) estimator
  - `analyze_feature_network()` - Network centrality metrics (degree, strength, betweenness)
  - `filter_by_network_connectivity()` - Remove disconnected features (likely artifacts)
  - `plot_feature_network()` - Visualize feature connectivity network

### Sparsity & Outlier Detection
- **`identify_sparsity_elbow.R`** (`identify_sparsity_elbow`) - Detect elbow point in richness-depth curve
- **`compute_depth_sparsity.R`** (`analyze_depth_sparsity`, `filter_depth_sparsity_outliers`) - Depth-sparsity outlier analysis
- **`plot_reads_vs_asvs.R`** - Total reads vs observed ASVs plot with outlier flagging

### Joint Filtering
- **`filter_features_joint.R`** (`filter_features_joint`) - Joint abundance-prevalence filtering with AND/OR logic

### Diversity Metrics
- **`calculate_cv.R`** (`calculate_feature_cv`) - Coefficient of variation per feature
- **`calculate_sparsity.R`** (`calculate_sparsity`) - Calculate sparsity metric

### QC & Metrics
- **`compute_qc_metrics.R`** (`compute_filtering_qc`) - Comprehensive QC including:
  - Sparsity changes (original vs filtered)
  - Retention rates (reads, features, samples)
  - Rank-abundance stability (top N overlap, Jaccard similarity, Spearman correlation)
  - Procrustes analysis for compositional similarity (via vegan)
  - Hill numbers / Effective Species Count diversity retention

### Visualization
- **`plot_qc_comparison.R`** (`plot_qc_comparison`) - Generates 7 comparison plots:
  1. Sample coverage distribution (log10, faceted)
  2. Feature abundance distribution (log10, faceted)
  3. Per-sample sparsity histogram
  4. Retention rates barplot
  5. Top features stacked barplot
  6. Heatmaps (original and filtered)
  7. Presence frequency comparison
- **`plot_presence_analysis.R`** (`plot_presence_analysis`) - Feature prevalence and sample richness histograms (with optional comparison mode)
- **`plot_coverage_histogram.R`** (`plot_coverage_histogram`) - Sample coverage histograms with threshold line
- **`plot_top_features_stacked.R`** (`plot_top_features_stacked`) - Top N relative abundance stacked bars
- **`plot_scree.R`** (`plot_scree`) - Scree/saturation diagnostic plots
- **`plot_depth_sparsity.R`** (`plot_depth_sparsity`) - Depth-sparsity relationship with outliers
- **`plot_sparsity_elbow.R`** (`plot_sparsity_elbow`) - Richness-depth curve with elbow marker

### Pipeline & Utilities
- **`run_filtering_pipeline.R`** (`run_filtering_pipeline`) - Complete filtering workflow with:
  - Accepts `data.frame`, `matrix`, `phyloseq`, or `TreeSummarizedExperiment` as input
  - Returns same class as input (preserving metadata)
  - Configurable methods per step (coverage, abundance)
  - File outputs (filtered table, plots, text report)
- **`load_feature_table.R`** (`load_feature_table`) - Auto-detects TSV/CSV format, header detection, feature column identification
- **`convert_to_relative.R`** (`convert_to_relative`) - Converts counts to relative abundances

### Data Class Conversion
- **`convert_data_classes.R`** - Native support for phyloseq and TreeSummarizedExperiment:
  - `from_phyloseq()` - Convert phyloseq object to data.frame
  - `to_phyloseq()` - Create phyloseq from data.frame (with optional taxonomy, tree, sample data)
  - `from_TSE()` - Convert TreeSummarizedExperiment/SingleCellExperiment to data.frame
  - `to_TSE()` - Create TreeSummarizedExperiment from data.frame
  - `convert_feature_table()` - Generic auto-detect conversion between any supported format

### Scree Analysis
- **`compute_scree.R`** (`compute_scree`) - Scree/saturation analysis for filtering threshold selection

### Key Design Patterns

1. **Feature Table Format**: First column contains feature IDs, remaining columns are sample counts. Functions preserve this structure (don't convert to row names).

2. **Method Dispatch**: Pipeline uses `match.arg()` for method selection:
   - Coverage: "none", "absolute", "mad", "iqr", "good", "chao"
   - Abundance: "none", "absolute", "relative", "relative_cutoff", "joint"

3. **Return Values**: Most functions return modified data.frames with attributes (e.g., `n_filtered_out`, `threshold`, `mode`).

4. **Optional Dependencies**: Uses `requireNamespace()` for optional packages (vegan, pheatmap, ComplexHeatmap, phyloseq, TreeSummarizedExperiment) with graceful degradation.

5. **QC Metrics Structure**: `compute_filtering_qc()` returns comprehensive list including sparsity changes, retention percentages, top N feature comparisons with overlap statistics, rank correlations, and Procrustes analysis.

## Testing Structure

Tests use testthat v3 with focused test files:
- `test-load_feature_table.R` - File loading and auto-detection
- `test-filter_coverage.R` - Sample coverage filtering
- `test-filter_features.R` - Feature abundance filtering
- `test-compute_qc_metrics.R` - QC metric calculations
- `test-helpers.R` - Utility function tests
- `test-estimate_coverage_metrics.R` - Coverage estimator tests
- `test-filter_by_singleton_ratio.R` - Singleton ratio filtering
- `test-filter_cross_talk.R` - Cross-talk filtering
- `test-filter_features_joint.R` - Joint filtering
- `test-compute_depth_sparsity.R` - Depth-sparsity analysis
- `test-identify_sparsity_elbow.R` - Elbow detection
- `test-compute_scree.R` - Scree analysis
- `test-filter_by_network.R` - Network-based filtering
- `test-convert_data_classes.R` - Class conversion functions

Test data is minimal synthetic data.frames created inline.

## Example Data

- `example_feature_table.tsv` - Synthetic dataset (50 features × 20 samples) for documentation and testing
