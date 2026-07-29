# featuretablefilter Tutorial

This tutorial walks you through installing the `featuretablefilter` package and running both basic and advanced filtering examples.

## Table of Contents

1. [Installation](#installation)
2. [Basic Example](#basic-example)
3. [Full Example](#full-example)
4. [Example Data](#example-data)

---

## Installation

### From Bioconductor (Recommended for Stable Release)

```r
# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install featuretablefilter from Bioconductor
BiocManager::install("featuretablefilter")
```

### From Source via devtools (Development Version)

For the latest development version, install using `devtools`:

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Install featuretablefilter from GitHub
devtools::install_github("VojtechBarton/featuretablefilter")

# Alternative: Install with dependencies if you encounter issues
devtools::install_github("VojtechBarton/featuretablefilter", dependencies = TRUE)
```

### Loading the Package

After installation, load the package:

```r
library(featuretablefilter)
```

---

## Basic Example

This example demonstrates a simple filtering workflow using the built-in example data.

```r
# Load the package
library(featuretablefilter)

# Load the built-in example data
data(example_feature_table)

# View the structure of the example data
str(example_feature_table)
head(example_feature_table)

# Run a basic filtering pipeline
result <- run_filtering_pipeline(
  input = example_feature_table,           # Use the example data directly
  output_dir = "basic_results",            # Directory for results
  prefix = "filtered",                     # Prefix for output files
  
  # Coverage filtering: remove low-coverage samples
  cov_filter_method = "mad",               # MAD-based method (robust to outliers)
  cov_threshold = 2,                       # Remove samples > 2 MAD below median
  
  # Abundance filtering: remove rare features
  abun_filter_method = "relative",         # Relative abundance
  abun_threshold = 0.001                   # Remove features < 0.1% relative abundance
)

# Access the results
filtered_table <- result$filtered_table
qc_metrics <- result$qc_metrics

# Compare original vs filtered dimensions
cat("Original table:", nrow(example_feature_table), "features x", 
    ncol(example_feature_table) - 1, "samples\n")
cat("Filtered table:", nrow(filtered_table), "features x", 
    ncol(filtered_table) - 1, "samples\n")

# View filtering summary
print(result$filtering_summary)
```

### What This Does

1. **Coverage Filtering**: Uses the Median Absolute Deviation (MAD) method to identify and remove samples with unusually low sequencing depth. This is a robust statistical approach that handles outliers well.

2. **Abundance Filtering**: Removes features (e.g., ASVs or OTUs) that have less than 0.1% relative abundance in any sample. This helps eliminate spurious sequences and sequencing errors.

3. **QC Metrics**: Automatically calculates comprehensive quality control metrics comparing before and after filtering.

### Expected Output

The pipeline creates:
- `filtered_filtered.tsv` - The filtered feature table
- `qc_report.txt` - Text report with filtering statistics
- Summary information in R including:
  - Number of samples/features removed
  - Retention rates
  - Sparsity changes
  - Diversity metrics

---

## Full Example

This comprehensive example demonstrates all major filtering methods, visualization options, and report generation.

```r
# Load required libraries
library(featuretablefilter)
library(ggplot2)          # For additional plotting (optional)

# ============================================================================
# Step 1: Load and Explore Example Data
# ============================================================================

data(example_feature_table)

# Basic exploration
cat("Table dimensions:", dim(example_feature_table), "\n")
cat("Feature IDs:", head(example_feature_table[[1]]), "\n")
cat("Sample IDs:", colnames(example_feature_table)[-1], "\n")

# Calculate basic statistics
total_reads <- sum(example_feature_table[, -1])
sample_depths <- colSums(example_feature_table[, -1])
feature_totals <- rowSums(example_feature_table[, -1])

cat("\nTotal reads:", total_reads, "\n")
cat("Sample depth range:", range(sample_depths), "\n")
cat("Feature count range:", range(feature_totals), "\n")

# ============================================================================
# Step 2: Coverage Filtering Exploration
# ============================================================================

# Method A: Absolute threshold
result_absolute <- run_filtering_pipeline(
  input = example_feature_table,
  cov_filter_method = "absolute",
  cov_threshold = 1000,      # Minimum 1000 reads per sample
  abun_filter_method = "none"
)

# Method B: MAD-based (recommended)
result_mad <- run_filtering_pipeline(
  input = example_feature_table,
  cov_filter_method = "mad",
  cov_threshold = 1.5,       # 1.5 MAD below median
  abun_filter_method = "none"
)

# Method C: Good's coverage estimator
result_good <- run_filtering_pipeline(
  input = example_feature_table,
  cov_filter_method = "good",
  cov_target_coverage = 0.95,  # Target 95% ecological coverage
  abun_filter_method = "none"
)

# Compare methods
cat("\n=== Coverage Filtering Comparison ===\n")
cat("Absolute (1000 reads):", result_absolute$qc_metrics$samples_retained, "samples retained\n")
cat("MAD (k=1.5):", result_mad$qc_metrics$samples_retained, "samples retained\n")
cat("Good's (95%):", result_good$qc_metrics$samples_retained, "samples retained\n")

# ============================================================================
# Step 3: Singleton Ratio Filtering
# ============================================================================

# Filter samples with excessive singletons (potential PCR artifacts)
result_singleton <- run_filtering_pipeline(
  input = example_feature_table,
  cov_filter_method = "none",
  singleton_max_ratio = 0.1,     # Max 10% singletons + doubletons
  singleton_count_type = "both", # Consider both singletons and doubletons
  abun_filter_method = "none"
)

# ============================================================================
# Step 4: Abundance Filtering
# ============================================================================

# Method A: Absolute count threshold
result_abs_abun <- run_filtering_pipeline(
  input = example_feature_table,
  cov_filter_method = "none",
  abun_filter_method = "absolute",
  abun_threshold = 10,           # Minimum 10 reads per feature
  abun_min_samples = 2           # Must appear in at least 2 samples
)

# Method B: Relative abundance
result_rel_abun <- run_filtering_pipeline(
  input = example_feature_table,
  cov_filter_method = "none",
  abun_filter_method = "relative",
  abun_threshold = 0.001         # 0.1% relative abundance
)

# Method C: Joint abundance-prevalence filtering
result_joint <- run_filtering_pipeline(
  input = example_feature_table,
  cov_filter_method = "none",
  abun_filter_method = "joint",
  abun_threshold = 0.001,        # 0.1% abundance
  abun_prevalence_threshold = 0.1, # Present in at least 10% of samples
  abun_logic = "AND"             # Feature must meet BOTH criteria
)

# ============================================================================
# Step 5: Cross-talk / Index-hopping Correction
# ============================================================================

# Identify and correct potential cross-contamination
result_crosstalk <- run_filtering_pipeline(
  input = example_feature_table,
  cov_filter_method = "none",
  abun_filter_method = "none",
  crosstalk_filter_method = "zero",      # Set suspected leakage to zero
  crosstalk_threshold = 0.001            # 0.1% of max abundance
)

# ============================================================================
# Step 6: Advanced Filtering Methods
# ============================================================================

# Sparsity elbow detection (diagnostic analysis)
elbow_result <- identify_sparsity_elbow(
  table = example_feature_table,
  method = "kneedle"
)

# Print filtering recommendation
print(elbow_result$recommendation)

# To apply elbow-based filtering in a pipeline:
result_elbow <- run_filtering_pipeline(
  input = example_feature_table,
  cov_filter_method = "mad",           # Use valid coverage method
  sparsity_elbow_detect = TRUE,        # Enable elbow detection
  apply_sparsity_elbow = FALSE,        # Set TRUE to apply elbow filtering
  sparsity_elbow_method = "kneedle",   # Elbow detection algorithm
  abun_filter_method = "none"
)

# Depth-sparsity outlier analysis
ds_analysis <- analyze_depth_sparsity(
  table = example_feature_table,
  metric = "sparsity",
  outlier_method = "mad",
  multiplier = 3,
  direction = "high_sparsity"
)

# View detected outliers
print(ds_analysis$outliers)

# Filter out depth-sparsity outliers
result_ds_filtered <- filter_depth_sparsity_outliers(
  table = example_feature_table,
  metric = "sparsity",
  outlier_method = "mad",
  multiplier = 3,
  direction = "high_sparsity"
)

# ============================================================================
# Step 7: Complete Pipeline with All Filters
# ============================================================================

# Comprehensive filtering combining multiple methods
result_full <- run_filtering_pipeline(
  input = example_feature_table,
  output_dir = "full_results",
  prefix = "filtered",
  
  # Sample-level filters
  cov_filter_method = "mad",
  cov_threshold = 2,
  singleton_max_ratio = 0.15,
  
  # Feature-level filters
  abun_filter_method = "joint",
  abun_threshold = 0.001,
  abun_prevalence_threshold = 0.1,
  abun_logic = "AND",
  
  # Cross-talk correction
  crosstalk_filter_method = "zero",
  crosstalk_threshold = 0.001,

  # Output options
  generate_plots = TRUE,
  generate_report = TRUE
)

# ============================================================================
# Step 8: Analyze Results
# ============================================================================

# View comprehensive QC metrics
qc <- result_full$qc_metrics

# Samples retained
cat("\n=== Sample Retention ===\n")
cat("Original table:", nrow(result_full$original_table), "features x",
    ncol(result_full$original_table) - 1, "samples\n")
cat("Filtered table:", nrow(result_full$filtered_table), "features x",
    ncol(result_full$filtered_table) - 1, "samples\n")
cat("Sample retention:", round(qc$sample_retention_percent, 1), "%\n")

# Features retained
cat("\n=== Feature Retention ===\n")
cat("Feature retention:", round(qc$feature_retention_percent, 1), "%\n")
cat("Read retention:", round(qc$read_retention_percent, 1), "%\n")

# Sparsity changes
cat("\n=== Sparsity Changes ===\n")
cat("Original sparsity:", round(qc$sparsity_original * 100, 1), "%\n")
cat("Filtered sparsity:", round(qc$sparsity_filtered * 100, 1), "%\n")
cat("Sparsity drop:", round(qc$sparsity_drop_percent, 1), "percentage points\n")

# Diversity comparison (Effective Number of Species / Hill numbers)
cat("\n=== Diversity Metrics ===\n")
cat("Shannon ENS - Original:", round(qc$shannon_ens_original, 2), "\n")
cat("Shannon ENS - Filtered:", round(qc$shannon_ens_filtered, 2), "\n")
cat("Shannon ENS retention:", round(qc$shannon_ens_retention_percent, 1), "%\n")
cat("Simpson ENS - Original:", round(qc$simpson_ens_original, 2), "\n")
cat("Simpson ENS - Filtered:", round(qc$simpson_ens_filtered, 2), "\n")
cat("Simpson ENS retention:", round(qc$simpson_ens_retention_percent, 1), "%\n")

# Rank-abundance stability
cat("\n=== Rank-Abundance Stability ===\n")
cat("Top 10 overlap:", qc$top_n_overlap_count, "features\n")
cat("Rank-abundance correlation:", round(qc$rank_abundance_correlation, 3), "\n")

# ============================================================================
# Step 9: Visualization
# ============================================================================

# Note: Visualization functions require ggplot2 (included via Imports)

# Coverage histogram with threshold line
# Wrap in print() to display in RStudio
print(plot_coverage_histogram(result_full$original_table)$plot)

# Before/after comparison plot (returns list of plots)
qc_plots <- plot_qc_comparison(
  original_table = result_full$original_table,
  filtered_table = result_full$filtered_table
)
# Print individual plots from the list
print(qc_plots$plots$coverage_distribution)
print(qc_plots$plots$sparsity_histogram)
print(qc_plots$plots$retention_rates)

# Top features stacked barplot (requires both tables)
print(plot_top_features_stacked(
  original_table = result_full$original_table,
  filtered_table = result_full$filtered_table,
  top_n = 10
))

# Optional: Sparsity elbow plot (if elbow detection was run)
if (!is.null(result_full$sparsity_elbow_result)) {
  print(plot_sparsity_elbow(result_full$sparsity_elbow_result))
}

# Optional: Scree plot (if scree analysis was run)
if (!is.null(result_full$scree_result)) {
  print(plot_scree(result_full$scree_result))
}

# ============================================================================
# Step 10: Export Results
# ============================================================================

# Save filtered table
write.table(
  result_full$filtered_table,
  file = "full_results/filtered_filtered.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Save QC metrics as CSV
write.csv(
  qc,
  file = "full_results/qc_metrics.csv",
  row.names = FALSE
)

# ============================================================================
# Working with phyloseq Objects
# ============================================================================

# If you have phyloseq installed, the package works directly with phyloseq objects:

# library(phyloseq)
# ps <- import_biom("your_data.biom")
# 
# result_ps <- run_filtering_pipeline(
#   input = ps,
#   cov_filter_method = "mad",
#   abun_filter_method = "relative",
#   abun_threshold = 0.001
# )
# 
# # Result is also a phyloseq object
# class(result_ps$filtered_table)  # "phyloseq"

# ============================================================================
# Working with TreeSummarizedExperiment
# ============================================================================

# If you have TreeSummarizedExperiment installed:

# library(TreeSummarizedExperiment)
# tse <- readRDS("your_data.rds")
# 
# result_tse <- run_filtering_pipeline(
#   input = tse,
#   cov_filter_method = "mad",
#   abun_filter_method = "relative",
#   abun_threshold = 0.001
# )
# 
# # Result preserves TreeSummarizedExperiment class
# class(result_tse$filtered_table)  # "TreeSummarizedExperiment"
```

---

## Example Data

The package includes `example_feature_table`, a synthetic dataset for testing and demonstration purposes.

### Dataset Specifications

| Property | Value |
|----------|-------|
| **Features** | 50 ASVs/OTUs |
| **Samples** | 20 samples |
| **Format** | data.frame |
| **Structure** | First column = feature IDs, remaining columns = sample counts |
| **Source** | Synthetic data generated for package demonstration |

### Loading the Data

```r
# Load the example data
data(example_feature_table)

# View structure
str(example_feature_table)

# View first few rows
head(example_feature_table)

# Check dimensions
dim(example_feature_table)

# Get sample names
sample_names <- colnames(example_feature_table)[-1]

# Get feature names
feature_names <- example_feature_table[[1]]
```

### Using Example Data in Your Analysis

```r
# Quick test of filtering parameters
test_result <- run_filtering_pipeline(
  input = example_feature_table,
  cov_filter_method = "mad",
  abun_filter_method = "relative",
  abun_threshold = 0.001
)

# Use as template for your own data
# Replace example_feature_table with your data file path:
# input = "path/to/your/feature_table.tsv"
```

### File Format Requirements

The example data follows the standard wide format expected by the package:

```
FeatureID    Sample1    Sample2    Sample3    ...
ASV001       156        89         234        ...
ASV002       0          45         12         ...
ASV003       78         0          156        ...
...
```

**Key points:**
- First column contains feature identifiers (no header requirement varies by import function)
- Remaining columns are sample counts
- Tab-separated (.tsv) or comma-separated (.csv) formats supported
- Zero counts should be explicit (not empty cells)

---

## Next Steps

After completing this tutorial, explore:

1. **Scree Analysis**: Systematically evaluate threshold effects
   ```r
   scree_result <- run_scree_analysis(example_feature_table, 
                                      cov_range = c(500, 5000),
                                      abun_range = c(0.0001, 0.01))
   ```

2. **Network-based Filtering**: Use mutual information networks
   ```r
   network_result <- apply_network_filtering(filtered_table)
   ```

3. **Interactive Dashboard**: Launch the Shiny dashboard for visual parameter exploration
   ```r
   runDashboard()
   ```

4. **Advanced Reporting**: Generate publication-ready PDF reports
   ```r
   generate_pdf_report(result_full, output_file = "filtering_report.pdf")
   ```

---

## Troubleshooting

### Common Issues

**Issue**: `devtools::install_github()` fails with dependency errors
```r
# Solution: Install dependencies explicitly
devtools::install_github("VojtechBarton/featuretablefilter", 
                         dependencies = c("Depends", "Imports", "Suggests"))
```

**Issue**: Functions not found after installation
```r
# Solution: Ensure package is loaded
library(featuretablefilter)

# Or use namespace qualification
result <- featuretablefilter::run_filtering_pipeline(...)
```

**Issue**: Missing optional visualization functions
```r
# Install required packages
install.packages(c("ggplot2", "pheatmap", "ComplexHeatmap", "vegan"))
```

---

## Additional Resources

- **Bioconductor Page**: https://bioconductor.org/packages/featuretablefilter
- **GitHub Repository**: https://github.com/VojtechBarton/featuretablefilter
- **Package README**: See `README.md` in the repository root
- **Support**: https://support.bioconductor.org/

---

*Last updated: July 2026*
