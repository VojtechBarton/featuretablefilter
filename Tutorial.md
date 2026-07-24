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

# Sparsity elbow detection
result_elbow <- run_filtering_pipeline(
  input = example_feature_table,
  cov_filter_method = "elbow",
  elbow_method = "kneedle",
  abun_filter_method = "none"
)

# Depth-sparsity outlier analysis
result_ds_outliers <- run_depth_sparsity_analysis(
  input = example_feature_table,
  metric = "sparsity",
  method = "mad"
)

# View detected outliers
print(result_ds_outliers$outliers)

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
  
  # Generate reports
  generate_text_report = TRUE,
  generate_markdown_report = TRUE
)

# ============================================================================
# Step 8: Analyze Results
# ============================================================================

# View comprehensive QC metrics
qc <- result_full$qc_metrics

# Samples retained
cat("\n=== Sample Retention ===\n")
cat("Original samples:", qc$original$n_samples, "\n")
cat("Filtered samples:", qc$filtered$n_samples, "\n")
cat("Retention rate:", round(qc$sample_retention_rate * 100, 1), "%\n")

# Features retained
cat("\n=== Feature Retention ===\n")
cat("Original features:", qc$original$n_features, "\n")
cat("Filtered features:", qc$filtered$n_features, "\n")
cat("Retention rate:", round(qc$feature_retention_rate * 100, 1), "%\n")

# Sparsity changes
cat("\n=== Sparsity Changes ===\n")
cat("Original sparsity:", round(qc$original$sparsity * 100, 1), "%\n")
cat("Filtered sparsity:", round(qc$filtered$sparsity * 100, 1), "%\n")

# Diversity comparison
cat("\n=== Diversity Metrics ===\n")
print(qc$diversity_comparison)

# ============================================================================
# Step 9: Visualization
# ============================================================================

# Note: Visualization functions require optional packages (pheatmap, ComplexHeatmap)
# Install with: install.packages(c("pheatmap", "ComplexHeatmap"))

# Coverage distribution plot
plot_coverage_distribution(example_feature_table)

# Before/after comparison
plot_qc_comparison(qc)

# Top features stacked barplot
plot_top_features_stacked(example_feature_table, top_n = 10)

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
