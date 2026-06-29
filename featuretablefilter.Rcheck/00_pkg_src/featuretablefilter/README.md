# featuretablefilter

An R package for filtering microbiome feature tables based on coverage and abundance criteria. The package provides functions for sample coverage filtering, feature abundance filtering, QC metrics computation, and visualization.

## Installation

```r
# Install from local source
devtools::install(".")
```

## Features

### Core Filtering Functions

- **`filter_by_coverage()`** - Removes low-coverage samples using absolute thresholds
- **`filter_features_by_abundance()`** - Filters features by absolute or relative abundance with `min_samples` support
- **`filter_by_relative_cutoff()`** - Relative threshold filtering based on min-coverage sample

### Coverage Estimation Functions

- **`estimate_mad_cutoff()`** - MAD-based threshold estimation (median +/- multiplier * MAD)
- **`estimate_iqr_cutoff()`** - IQR/Tukey-based threshold estimation
- **`estimate_good_coverage()`** - Good's coverage estimator for ecological completeness
- **`estimate_chao_coverage()`** - Chao's coverage estimator (conservative completeness)
- **`filter_by_coverage_estimator()`** - Filter samples by ecological coverage threshold

### QC & Metrics

- **`compute_filtering_qc()`** - Computes sparsity, retention rates, rank-abundance stability (top N overlap, Jaccard similarity, Spearman correlation), and Procrustes analysis
- **`calculate_sparsity()`** - Calculate sparsity metric
- **`calculate_cv()`** - Calculate coefficient of variation

### Visualization

- **`plot_qc_comparison()`** - Generates comprehensive QC comparison plots:
  - Sample coverage distribution
  - Feature abundance distribution
  - Sample sparsity histogram (per-sample view)
  - Retention rates
  - Heatmaps of top variable features
  - Stacked barplots of top features
- **`plot_presence_analysis()`** - Feature prevalence and sample richness histograms with comparison capability
- **`plot_coverage_histogram()`** - Sample coverage histograms
- **`plot_top_features_stacked()`** - Top N relative abundance stacked bars

### Pipeline & Utilities

- **`run_filtering_pipeline()`** - Orchestrates complete filtering workflow with configurable methods per step, saves outputs, generates text reports
- **`load_feature_table()`** - Auto-detects TSV/CSV format, handles header detection, feature column identification
- **`convert_to_relative()`** - Converts counts to relative abundances

## Quick Start

### Simple Filtering Pipeline

```r
library(featuretablefilter)

# Run a complete filtering pipeline
result <- run_filtering_pipeline(
  input_file = "feature-table.tsv",
  output_dir = "results",
  prefix = "analysis1",
  cov_filter_method = "absolute",
  cov_threshold = 1000,
  abun_filter_method = "absolute",
  abun_threshold = 5
)
```

### Individual Functions

```r
# Load feature table
table <- load_feature_table("feature-table.tsv")

# Filter by coverage (remove samples with < 1000 reads)
filtered <- filter_by_coverage(table, min_reads = 1000)

# Filter features by abundance (remove features with < 5 total reads)
filtered <- filter_features_by_abundance(
  filtered,
  threshold = 5,
  mode = "absolute",
  min_samples = 1
)

# Compute QC metrics
qc <- compute_filtering_qc(original_table, filtered_table)

# Generate comparison plots
plot_qc_comparison(original_table, filtered_table, plot_dir = "plots")

# Presence frequency analysis with comparison
plot_presence_analysis(
  original_table,
  table_filtered = filtered_table,
  save_dir = "plots",
  prefix = "presence"
)

# Filter by ecological completeness (Good's coverage >= 95%)
filtered <- filter_by_coverage_estimator(
  table,
  method = "good",
  target_coverage = 0.95
)
cat(sprintf("Mean coverage before: %.2f%%\n", filtered$mean_coverage_before * 100))
cat(sprintf("Mean coverage after: %.2f%%\n", filtered$mean_coverage_after * 100))
```

## Filtering Methods

### Coverage Filtering (`cov_filter_method`)

| Method | Description |
|--------|-------------|
| `none` | No coverage filtering |
| `absolute` | Remove samples below fixed read count threshold |
| `mad` | Use MAD-based threshold (median +/- multiplier * MAD) |
| `iqr` | Use IQR/Tukey-based threshold (fences) |
| `good` | Use Good's coverage estimator (ecological completeness, default 95%) |
| `chao` | Use Chao's coverage estimator (conservative, default 90%) |

### Abundance Filtering (`abun_filter_method`)

| Method | Description |
|--------|-------------|
| `none` | No abundance filtering |
| `absolute` | Remove features below fixed count threshold |
| `relative` | Remove features below relative abundance proportion |
| `relative_cutoff` | Relative threshold based on min-coverage sample |

## Output Files

When running the pipeline, the following files are generated:

- `<prefix>_table.tsv` - Filtered feature table
- `<prefix>_table.tsv` - Filtered feature table
- `<prefix>_sample_coverage_distribution.png` - Coverage distribution histograms
- `<prefix>_feature_abundance_distribution.png` - Feature abundance histograms
- `<prefix>_sparsity_comparison.png` - Per-sample sparsity histograms
- `<prefix>_retention_rates.png` - Retention rate barplot
- `<prefix>_heatmap_original.png` - Top variable features heatmap (original)
- `<prefix>_heatmap_filtered.png` - Top variable features heatmap (filtered)
- `<prefix>_top_features_stacked.png` - Top features stacked barplot
- `<prefix>_presence_frequency_comparison.png` - Prevalence/richness comparison
- `<prefix>_report.txt` - Text summary report

## License

GPL-3
