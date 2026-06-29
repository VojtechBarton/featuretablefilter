# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

**featuretablefilter** - An R package for filtering microbiome feature tables based on coverage and abundance criteria. The package provides functions for sample coverage filtering, feature abundance filtering, QC metrics computation, and visualization.

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
```

## Architecture

The package follows a modular structure with distinct functional categories:

### Core Filtering Functions
- **`filter_by_coverage.R`** - Removes low-coverage samples using absolute thresholds
- **`filter_features.R`** - Filters features by absolute or relative abundance with `min_samples` support
- **`filter_relative_cutoff.R`** - Relative cutoff filtering based on min-coverage sample

### Coverage Estimation Functions
- **`estimate_mad_cutoff.R`** - MAD-based threshold estimation (median ± multiplier × MAD)
- **`estimate_iqr_cutoff.R`** - IQR/Tukey-based threshold estimation

### QC & Metrics
- **`compute_qc_metrics.R`** (`compute_filtering_qc`) - Computes sparsity, retention rates, rank-abundance stability (top N overlap, Jaccard similarity, Spearman correlation), and Procrustes analysis (via vegan)
- **`calculate_sparsity.R`**, **`calculate_cv.R`** - Helper metrics

### Visualization
- **`plot_qc_comparison.R`** (`plot_qc_comparison`) - Generates 7 plots: coverage distribution, feature abundance distribution, sparsity comparison, retention rates, heatmaps, stacked barplots
- **`plot_presence_analysis.R`** - Feature prevalence histograms
- **`plot_coverage_histogram.R`** - Sample coverage histograms
- **`plot_top_features_stacked.R`** - Top N relative abundance stacked bars

### Pipeline & Utilities
- **`run_filtering_pipeline.R`** (`run_filtering_pipeline`) - Orchestrates complete filtering workflow with configurable methods per step, saves outputs, generates text reports
- **`load_feature_table.R`** (`load_feature_table`) - Auto-detects TSV/CSV format, handles header detection, feature column identification
- **`convert_to_relative.R`** - Converts counts to relative abundances

### Key Design Patterns

1. **Feature Table Format**: First column contains feature IDs, remaining columns are sample counts. Functions preserve this structure (don't convert to row names).

2. **Method Dispatch**: Pipeline uses `match.arg()` for method selection ("none", "absolute", "mad", "iqr" for coverage; "none", "absolute", "relative", "relative_cutoff" for abundance).

3. **Return Values**: Most functions return modified data.frames with attributes (e.g., `n_filtered_out`, `threshold`, `mode`).

4. **Optional Dependencies**: Uses `requireNamespace()` for optional packages (vegan, pheatmap, ComplexHeatmap) with graceful degradation.

5. **QC Metrics Structure**: `compute_filtering_qc()` returns comprehensive list including sparsity changes, retention percentages, top N feature comparisons with overlap statistics and rank correlations.

## Testing Structure

Tests use testthat v3 with focused test files:
- `test-load_feature_table.R` - File loading and auto-detection
- `test-filter_coverage.R` - Sample coverage filtering
- `test-filter_features.R` - Feature abundance filtering
- `test-compute_qc_metrics.R` - QC metric calculations
- `test-helpers.R` - Utility function tests

Test data is minimal synthetic data.frames created inline.
