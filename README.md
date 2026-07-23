# featuretablefilter

[![Bioconductor-time](http://bioconductor.org/shields/years-in-bioc/featuretablefilter.svg)](https://bioconductor.org/packages/release/bioc/html/featuretablefilter.html)
[![Bioconductor Downloads](http://bioconductor.org/shields/downloads/release/featuretablefilter.svg)](https://bioconductor.org/packages/stats/bioc/featuretablefilter/)

An R package for filtering microbiome feature tables based on coverage, abundance, and network connectivity criteria. Provides comprehensive quality control metrics and visualizations for filtering decisions. **Native support for `phyloseq` and `TreeSummarizedExperiment` objects.**

## Installation

### From Bioconductor

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("featuretablefilter")
```

### From Source (Development)

```r
# Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Install from GitHub
devtools::install_github("VojtechBarton/featuretablefilter")
```

## Quick Start

```r
library(featuretablefilter)

# Run a complete filtering pipeline
result <- run_filtering_pipeline(
  input = "path/to/feature_table.tsv",
  output_dir = "results",
  prefix = "filtered",
  cov_filter_method = "mad",       # MAD-based coverage filtering
  cov_threshold = 2,               # 2 MAD below median
  abun_filter_method = "relative", # Relative abundance filtering
  abun_threshold = 0.001           # 0.1% relative abundance
)

# Access results
filtered_table <- result$filtered_table
qc_metrics <- result$qc_metrics
```

## Features

### Core Filtering

#### Coverage Filtering (Sample-level)

Remove low-coverage samples that may introduce noise or bias. Available methods:

| Method | Description | When to Use |
|--------|-------------|-------------|
| `absolute` | Fixed minimum read count threshold | When you have a known minimum depth requirement |
| `mad` | MAD-based threshold (median ± k × MAD) | Default choice; robust to outliers |
| `iqr` | Tukey's IQR method (Q1 - k × IQR) | Alternative robust method; similar to MAD |
| `good` | Good's coverage estimator | When ecological completeness is priority |
| `chao` | Chao's coverage estimator | More conservative than Good's |

**Parameters:**
- `cov_threshold`: Fixed value for `absolute`, multiplier for `mad`/`iqr` (default ~1.5-2)
- `cov_floor`: Minimum coverage floor when using MAD/IQR
- `cov_target_coverage`: Target coverage (0-1) for Good's/Chao's methods

#### Singleton Ratio Filtering (Sample-level)

Detect and remove samples with excessive singletons/doubletons, which often indicate poor sequencing quality or PCR artifacts.

**Parameters:**
- `singleton_max_ratio`: Maximum allowed ratio of (singletons + doubletons) / total reads (default 0.1 = 10%)
- `singleton_count_type`: `"singleton"`, `"doubleton"`, or `"both"`

#### Abundance Filtering (Feature-level)

Filter features based on their abundance across samples. Available methods:

| Method | Description | When to Use |
|--------|-------------|-------------|
| `absolute` | Minimum raw read count per feature | When working with count data and have fixed cutoff |
| `relative` | Minimum relative abundance proportion | For compositional comparisons across samples |
| `relative_cutoff` | Relative threshold based on min-coverage sample | Adaptive; accounts for varying library sizes |
| `joint` | Combined abundance AND/OR prevalence criteria | Stringent filtering; requires both abundance and presence |

**Parameters:**
- `abun_threshold`: Minimum abundance (count or proportion)
- `abun_min_samples`: Minimum samples where feature must exceed threshold
- `abun_logic`: `"AND"` or `"OR"` for joint filtering
- `abun_prevalence_threshold`: Proportion of samples for joint filtering

#### Cross-talk / Index-hopping Filtering (Feature-level)

Identify and correct cross-contamination artifacts from Illumina patterned flow cells.

**Methods:**

| Method | Description |
|--------|-------------|
| `zero` | Set suspected leakage reads to zero (default) |
| `remove_feature` | Remove entire feature if any leakage detected |
| `flag` | Flag leakage but don't modify data |

**Parameters:**
- `crosstalk_threshold`: Maximum relative abundance for leakage detection (default 0.001 = 0.1% of max)
- `crosstalk_min_abs_cutoff`: Optional absolute count override

### Advanced Filtering

#### Sparsity Elbow Detection

Detects the elbow point in richness-depth curves to recommend optimal coverage cutoffs. Can be used diagnostically or applied as a filter.

**Methods:**
- `kneedle`: Kneedle algorithm for elbow detection
- `max_derivative`: Maximum derivative method
- `second_derivative`: Second derivative zero-crossing

#### Depth-Sparsity Outlier Analysis

Identifies samples that are outliers in the depth-sparsity relationship, which may indicate technical artifacts or biological anomalies.

**Metrics:**
- `"sparsity"`: Proportion of zeros per sample
- `"richness"`: Number of observed features per sample

**Outlier Detection:**
- `mad`: MAD-based outlier detection
- `iqr`: IQR-based outlier detection

#### Network-based Filtering

Uses mutual information networks to identify spurious features based on connectivity patterns.

#### Scree Analysis

Systematically evaluates filtering threshold effects across a range of values, helping to choose optimal parameters.

### Quality Control

- Comprehensive QC metrics including sparsity changes, retention rates, diversity retention
- Hill numbers (Effective Number of Species) for Shannon and Simpson diversity
- Procrustes analysis for compositional similarity
- Rank-abundance stability metrics (Jaccard similarity, Spearman correlation)

### Visualization

- Coverage and abundance distribution plots
- Sparsity elbow detection visualization
- Depth-sparsity outlier plots
- Top features stacked barplots
- Heatmaps and presence frequency plots
- Complete QC comparison reports

### Data Format Support

- **Input**: `data.frame`, `phyloseq`, `TreeSummarizedExperiment`, TSV/CSV files
- **Output**: Same class as input (preserves metadata and tree structures)

## Example Data

The package includes `example_feature_table`, a synthetic dataset (50 features × 20 samples) for testing and demonstration:

```r
data(example_feature_table)
str(example_feature_table)
```

## Reporting

Generate comprehensive reports in multiple formats:

- **Text report** - Plain text summary with all filtering steps and QC metrics
- **Markdown report** - Formatted report with tables, suitable for version control
- **PDF report** - Professional publication-ready document (requires LaTeX)

## Interactive Dashboard

Launch an interactive Shiny dashboard for exploring filtering parameters:

```r
runDashboard()
```

## Documentation

Full documentation is available:

- **Bioconductor**: https://bioconductor.org/packages/featuretablefilter
- **Reference manual**: `vignette("featuretablefilter")`
- **Tutorial**: `vignette("tutorial")`

## Citation

To cite this package in publications:

```r
citation("featuretablefilter")
```

## License

This package is licensed under the GNU GPL 3.0.

## Contributing

Please read our [Contributing Guide](CONTRIBUTING.md) before submitting pull requests.

## Support

For questions and issues:

- **Bioconductor Support Forum**: https://support.bioconductor.org/
- **GitHub Issues**: https://github.com/VojtechBarton/featuretablefilter/issues

---

*Built for the microbiome research community.*
