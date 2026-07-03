# featuretablefilter

An R package for filtering microbiome feature tables. Provides comprehensive tools for coverage-based filtering, abundance filtering, network-based filtering, QC metrics computation, and visualization. **Native support for phyloseq and TreeSummarizedExperiment objects.**

## Installation

```r
# Install from source
devtools::install(".")

# Or install from Bioconductor (when available)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("featuretablefilter")
```

## Features

### Core Filtering Functions

- **`filter_by_coverage()`** - Removes low-coverage samples using absolute thresholds
- **`filter_features_by_abundance()`** - Filters features by absolute or relative abundance with `min_samples` support
- **`filter_by_relative_cutoff()`** - Relative threshold filtering based on min-coverage sample

### Coverage Estimation Functions

- **`estimate_mad_cutoff()`** - MAD-based threshold estimation (median ± multiplier × MAD)
- **`estimate_iqr_cutoff()`** - IQR/Tukey-based threshold estimation
- **`estimate_good_coverage()`** - Good's coverage estimator for ecological completeness
- **`estimate_chao_coverage()`** - Chao's coverage estimator (conservative completeness)
- **`filter_by_coverage_estimator()`** - Filter samples by ecological coverage threshold

### Network-Based Filtering (New!)

Advanced filtering using information theory and graph theory to identify spurious features:

- **`compute_mutual_information()`** - Computes pairwise mutual information using k-nearest neighbor (KSG) estimator
- **`analyze_feature_network()`** - Calculates network centrality metrics (degree, strength, betweenness)
- **`filter_by_network_connectivity()`** - Removes features with no significant connections (likely artifacts)
- **`plot_feature_network()`** - Visualizes feature connectivity network

### Sparsity Elbow Detection (New!)

Identifies critical sequencing depth where ASV discovery crashes:

- **`identify_sparsity_elbow()`** - Detects elbow point using kneedle, max_derivative, or second_derivative methods
- **`plot_sparsity_elbow()`** - Visualizes richness-depth relationship with elbow marker
- **`analyze_depth_sparsity()`** - Full depth-sparsity analysis with MAD/IQR outlier detection
- **`plot_depth_sparsity()`** - Visualizes depth-sparsity relationship
- **`filter_depth_sparsity_outliers()`** - Convenience function to remove outlier samples
- **`plot_reads_vs_asvs()`** - Simple Total Reads vs Observed ASVs plot with MAD outlier flagging

### Diversity Retention Metrics (New!)

Effective Number of Species (Hill numbers) to assess how filtering affects biological diversity:

- **`calc_shannon_ens()`** - Shannon diversity converted to effective species count (Hill q=1)
- **`calc_simpson_ens()`** - Simpson diversity as effective species count (Hill q=2)
- **`calc_hill_numbers()`** - Full diversity profile across multiple q values
- **`compute_filtering_qc()`** - Now includes ENS retention percentages alongside traditional metrics

### Phyloseq & TreeSummarizedExperiment Support (New!)

Native conversion functions for standard microbiome data classes:

- **`from_phyloseq()`** - Convert phyloseq object to data.frame
- **`to_phyloseq()`** - Create phyloseq from data.frame (with optional taxonomy, tree, sample data)
- **`from_TSE()`** - Convert TreeSummarizedExperiment/SingleCellExperiment to data.frame
- **`to_TSE()`** - Create TreeSummarizedExperiment from data.frame
- **`convert_feature_table()`** - Generic auto-detect conversion between any supported format

### QC & Metrics

- **`compute_filtering_qc()`** - Comprehensive QC: sparsity, retention rates, rank-abundance stability, Procrustes analysis, **and ENS diversity retention**
- **`calculate_sparsity()`** - Calculate sparsity metric
- **`calculate_cv()`** - Calculate coefficient of variation

### Visualization

- **`plot_qc_comparison()`** - Comprehensive QC comparison plots (coverage, abundance, sparsity, retention, heatmaps, stacked bars)
- **`plot_presence_analysis()`** - Feature prevalence and sample richness histograms
- **`plot_coverage_histogram()`** - Sample coverage histograms
- **`plot_top_features_stacked()`** - Top N relative abundance stacked bars
- **`plot_scree()`** - Scree/saturation diagnostic plots for filtering thresholds
- **`plot_depth_sparsity()`** - Depth-sparsity relationship with outliers
- **`plot_sparsity_elbow()`** - Richness-depth curve with elbow detection

### Pipeline & Utilities

- **`run_filtering_pipeline()`** - Complete filtering workflow with configurable methods, file outputs, and reports
- **`load_feature_table()`** - Auto-detects TSV/CSV format, header detection, feature column identification
- **`convert_to_relative()`** - Converts counts to relative abundances

### Interactive Shiny Dashboard (New!)

**`runDashboard()`** - Launch an interactive web-based GUI for exploratory filtering without writing code.

## Quick Start

### Simple Filtering Pipeline

```r
library(featuretablefilter)

# Run a complete filtering pipeline
result <- run_filtering_pipeline(
  input = "example_feature_table.tsv",
  output_dir = "results",
  prefix = "analysis1",
  cov_filter_method = "absolute",
  cov_threshold = 1000,
  abun_filter_method = "absolute",
  abun_threshold = 5
)
```

### Working with phyloseq Objects

```r
library(phyloseq)
library(featuretablefilter)

# Load your phyloseq object
ps <- readRDS("my_phyloseq.rds")

# Convert to data.frame for filtering (optionally include taxonomy)
df <- from_phyloseq(ps, include_taxa = TRUE)

# Apply network-based filtering to remove spurious features
filtered_df <- filter_by_network_connectivity(df, similarity_type = "mi")

# Convert back to phyloseq (preserves tree if present)
ps_filtered <- to_phyloseq(filtered_df, phy_tree = phy_tree(ps))
```

### Working with TreeSummarizedExperiment

```r
library(TreeSummarizedExperiment)
library(featuretablefilter)

# Load TSE object
tse <- readRDS("my_tse.rds")

# Convert to data.frame with rowData and colData
result <- from_TSE(tse, add_row_data = TRUE, add_col_data = TRUE)
df <- result$table
sample_meta <- result$sample_data

# Filter and convert back
filtered_df <- filter_by_coverage(df, min_reads = 1000)
tse_filtered <- to_TSE(filtered_df, rowData = rowData(tse), colData = colData(tse))
```

### Network-Based Filtering

```r
# Remove features with no network connections (likely artifacts)
cleaned <- filter_by_network_connectivity(
  my_table,
  similarity_type = "mi",    # or "cor" for faster Pearson correlation
  min_degree = 1,            # minimum connections required
  verbose = TRUE
)

# Get detailed network analysis
mi_matrix <- compute_mutual_information(my_table)
network <- analyze_feature_network(mi_matrix, method = "mi")

# View disconnected features
disconnected <- names(network$degree_centrality[network$degree_centrality == 0])

# Visualize the network
plot <- plot_feature_network(network, rownames(my_table))
print(plot)
```

### Sparsity Elbow Detection

```r
# Identify critical depth where ASV discovery crashes
elbow_result <- identify_sparsity_elbow(
  my_table,
  method = "kneedle",        # or "max_derivative", "second_derivative"
  smooth_window = 5
)

cat(elbow_result$recommendation)
# "Consider filtering out 12 samples (24.0%) below depth 850 reads..."

# Visualize
plot_sparsity_elbow(elbow_result)

# Or use simple reads-vs-ASVs plot
read_asv_result <- plot_reads_vs_asvs(my_table, mad_multiplier = 3)
print(read_asv_result$plot)
outlier_samples <- read_asv_result$outliers$sample_name
```

### Diversity Retention Analysis

```r
# Compare original and filtered tables
qc <- compute_filtering_qc(original_table, filtered_table)

# Check how much biological diversity was retained
cat(sprintf("Shannon ENS retention: %.1f%%\n", qc$shannon_ens_retention_percent))
cat(sprintf("Simpson ENS retention: %.1f%%\n", qc$simpson_ens_retention_percent))

# Calculate Hill numbers directly
hill_profile <- calc_hill_numbers(as.matrix(table[, -1]), q = c(0, 0.5, 1, 2, 3))
```

### Interactive Shiny Dashboard

For those who prefer a graphical interface or want to explore different filtering strategies interactively:

```r
library(featuretablefilter)

# Launch the dashboard (opens in default browser)
runDashboard()

# Launch on specific port
runDashboard(port = 3838)

# Launch without auto-opening browser
runDashboard(launch.browser = FALSE)
```

**Features:**
- **Upload** - Load TSV/CSV feature tables directly in the browser
- **Interactive Filtering** - Configure coverage, singleton ratio, cross-talk, and abundance filtering with real-time preview
- **Visual Feedback** - See immediate updates to coverage distributions, abundance histograms, sparsity plots, and scree analysis
- **Retention Statistics** - Compare original vs. filtered data with retention rate visualizations
- **Export Options** - Download filtered tables (TSV), text reports, and visualization packages (ZIP)
- **Reproducible Code** - Generate R code for your filtering choices to reproduce in scripts

**Dashboard Tabs:**
1. **Overview** - Summary statistics and retention plot
2. **Coverage Distribution** - Sample coverage histograms
3. **Feature Abundance** - Feature abundance distribution
4. **Sparsity** - Per-sample sparsity analysis
5. **Scree Analysis** - Diagnostic plots for threshold selection
6. **Results** - Filtered table preview and download options
7. **R Code** - Reproducible code generation

## Individual Filtering Functions

```r
# Load feature table
table <- load_feature_table("example_feature_table.tsv")

# Filter by coverage (remove samples with < 1000 reads)
filtered <- filter_by_coverage(table, min_reads = 1000)

# Filter features by abundance (remove features with < 5 total reads)
filtered <- filter_features_by_abundance(
  filtered,
  threshold = 5,
  mode = "absolute",
  min_samples = 1
)

# Compute QC metrics (includes diversity retention)
qc <- compute_filtering_qc(original_table, filtered_table)

# Generate comparison plots
plot_qc_comparison(original_table, filtered_table, plot_dir = "plots")

# Presence frequency analysis
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
| `mad` | Use MAD-based threshold (median ± multiplier × MAD) |
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
| `joint` | Joint abundance-prevalence filtering with AND/OR logic |

## Output Files (Pipeline)

When running `run_filtering_pipeline()`, the following files are generated:

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

## Supported Data Classes

The package natively supports conversion between:

| Format | Description |
|--------|-------------|
| `data.frame` | Standard R data frame (feature ID in first column) |
| `phyloseq` | Legacy microbiome standard (OTU table, taxonomy, tree, sample data) |
| `TreeSummarizedExperiment` | Modern Bioconductor standard (assays, rowData, colData, reducedDims) |
| `SingleCellExperiment` | Single-cell variant of SummarizedExperiment |

## License

GPL-3
