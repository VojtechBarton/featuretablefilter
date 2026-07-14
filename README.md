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

- **Coverage filtering** - Remove low-coverage samples using absolute, MAD, IQR, or ecological coverage estimators (Good's, Chao's)
- **Abundance filtering** - Filter features by absolute count, relative abundance, or joint abundance-prevalence criteria
- **Singleton ratio filtering** - Detect and remove samples with excessive PCR artifacts
- **Cross-talk filtering** - Identify and correct index hopping contamination

### Advanced Filtering

- **Network-based filtering** - Use mutual information and network centrality to identify spurious features
- **Depth-sparsity analysis** - Detect outlier samples using richness-depth relationships
- **Scree analysis** - Systematically evaluate filtering threshold effects

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

This package is licensed under the Artistic License 2.0.

## Contributing

Please read our [Contributing Guide](CONTRIBUTING.md) before submitting pull requests.

## Support

For questions and issues:

- **Bioconductor Support Forum**: https://support.bioconductor.org/
- **GitHub Issues**: https://github.com/VojtechBarton/featuretablefilter/issues

---

*Built for the microbiome research community.*
