# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **`plot_presence_analysis()` comparison mode**: New `table_filtered` parameter enables side-by-side comparison of original vs filtered tables with:
  - Faceted histograms showing distributions for both datasets
  - Shared x-axis scales for fair visual comparison
  - Mean and median vertical lines with clear legend (dashed = mean, dot-dash = median)
  - Returns additional statistics: `mean_feature_prevalence_filtered`, `mean_sample_richness_filtered`, min/max values for filtered data

- **Per-sample sparsity histogram**: `plot_qc_comparison()` now shows sparsity distribution across samples instead of single aggregate values:
  - Histograms with faceting by dataset (Original/Filtered)
  - Mean sparsity vertical lines for each dataset
  - Better visualization of sparsity variation across samples

- **Enhanced heatmap documentation**: Heatmaps now include:
  - Clear title indicating rows = features, columns = samples
  - Dynamic row font sizing based on number of features
  - Larger output dimensions (1.5x) for better label visibility
  - Improved Z-score legend with proper labeling

- **`plot_presence_analysis()` histogram bins**: Changed from barplots to proper histograms with configurable bin counts (`prev_bins`, `rich_bins` parameters)

### Changed

- **Fixed feature prevalence / sample richness calculation**: Corrected swapped `rowSums`/`colSums` in `plot_presence_analysis()`:
  - `feature_prevalence` now correctly counts samples per feature (using `rowSums`)
  - `sample_richness` now correctly counts features per sample (using `colSums`)

- **Updated `run_filtering_pipeline()`**: Now uses comparison mode for presence analysis, generating single `<prefix>_presence_frequency_comparison.png` instead of separate plots

- **Improved histogram styling**: Thicker bin borders (linewidth 0.5) for better visibility

### Fixed

- Bug where y-axis showed incorrect counts (>3000 for 468 samples) due to swapped calculations
- Missing legend entries for mean/median statistics in comparison plots
- Duplicate vertical lines appearing in facet panels

## [0.1.0] - 2026-06-22

### Added

- Initial release with core filtering functions
- Coverage estimation methods (MAD, IQR)
- QC metrics computation
- Visualization functions
- Complete filtering pipeline with report generation
