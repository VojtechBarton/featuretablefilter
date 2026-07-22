# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

**featuretablefilter** is an R package for filtering microbiome feature tables by coverage, abundance, and network connectivity. It supports `data.frame`, `matrix`, `phyloseq`, and `TreeSummarizedExperiment` inputs and returns the same class as the input.

The package is organized as a standard R package with roxygen2 documentation, testthat v3 tests, Bioconductor metadata, and a new Shiny dashboard in `inst/shiny-dashboard/`.

## Development Commands

```bash
# Run all tests
Rscript -e "testthat::test_dir('tests/testthat')"

# Run a single test file
Rscript -e "testthat::test_file('tests/testthat/test-compute_qc_metrics.R')"

# Run tests matching a name filter
Rscript -e "testthat::test_file('tests/testthat/test-compute_qc_metrics.R', filter = 'sparsity')"

# Load package in development mode
Rscript -e "devtools::load_all('.')"

# Generate roxygen2 documentation
Rscript -e "roxygen2::roxygenise()"

# Build and install locally
R CMD INSTALL .

# Build source tarball
R CMD build .

# Run R package check (closest thing to linting for R packages)
R CMD check .

# Run the Shiny dashboard locally
Rscript -e "featuretablefilter::runDashboard()"
```

## High-Level Architecture

### Data Model

All functions assume a feature table in **wide format**: the first column contains feature IDs (e.g., ASVs, OTUs), and remaining columns contain per-sample counts. Functions generally preserve this structure rather than converting the first column to row names.

### Filtering Layers

Filtering is organized into independent, composable layers:

1. **Sample-level filtering** – remove low-coverage or low-quality samples.
   - Absolute thresholds, MAD/IQR estimators, Good’s/Chao’s ecological coverage.
   - Singleton-ratio filtering for PCR artifacts.
   - Cross-talk / index-hopping correction.
   - Depth-sparsity and sparsity-elbow outlier detection.

2. **Feature-level filtering** – remove rare or spurious features.
   - Absolute or relative abundance thresholds.
   - Relative cutoff tied to the minimum-coverage sample.
   - Joint abundance-prevalence filtering with AND/OR logic.
   - Network connectivity filtering based on mutual information.

3. **QC and visualization** – compare tables before and after filtering.
   - Retention rates, sparsity changes, diversity retention, rank-abundance stability.
   - Procrustes analysis via `vegan`.
   - Plotting functions for coverage, abundance, sparsity, scree, and QC comparison.

### Pipeline Internals

`run_filtering_pipeline()` is the main entry point. It delegates to internal pipeline modules in `R/`:

- `pipeline_steps.R` – wrappers such as `.apply_coverage_filter()`, `.apply_abundance_filter()`, `.apply_crosstalk_filter()`, `.apply_sparsity_elbow_filter()`, `.apply_depth_sparsity_filter()`. These standardize arguments and call the public filtering functions.
- `pipeline_io.R` – loading and saving inputs/outputs, preserving class.
- `pipeline_analysis.R` – orchestrates QC metrics, presence stats, scree, and depth-sparsity analysis.
- `pipeline_reporting.R` – generates text, Markdown, and PDF reports.

The pipeline returns a list with `original_table`, `filtered_table`, `qc_metrics`, `presence_stats`, `filtering_summary`, and optional analysis results.

### Class Conversion

`convert_data_classes.R` handles conversions to/from `phyloseq` and `TreeSummarizedExperiment`. The pipeline uses `convert_feature_table()` to coerce inputs to `data.frame` for processing, then converts back to the original class on output, preserving metadata and tree structures when possible.

### Key Design Patterns

- **Method dispatch** via `match.arg()`. Common method sets:
  - Coverage: `"none"`, `"absolute"`, `"mad"`, `"iqr"`, `"good"`, `"chao"`
  - Abundance: `"none"`, `"absolute"`, `"relative"`, `"relative_cutoff"`, `"joint"`
  - Cross-talk: `"none"`, `"zero"`, `"remove_feature"`, `"flag"`
- **Return values**: most functions return a modified `data.frame` with attributes such as `n_filtered_out`, `threshold`, and `mode`.
- **Optional dependencies**: use `requireNamespace()` for `vegan`, `pheatmap`, `ComplexHeatmap`, `phyloseq`, `TreeSummarizedExperiment`, `shiny`, and `DT` with graceful degradation.
- **Reports**: the pipeline can write a text report, a Markdown report, and a PDF report (PDF requires LaTeX).

## Shiny Dashboard

An interactive dashboard is shipped under `inst/shiny-dashboard/` and launched with `runDashboard()`.

- `app.R` – entry point.
- `ui.R` and `server.R` – sidebar parameters and reactive filtering/visualization.
- `www/style.css` – custom styling.

The dashboard supports file upload, real-time parameter exploration, retention-rate previews, dynamic plots, and export of the filtered table, QC report, plots, and reproducible R code. Optional Suggested dependencies: `shiny`, `DT`.

## Testing

Tests live in `tests/testthat/` and use testthat edition 3. The `helper.R` file loads `testthat` and the installed package. Synthetic `data.frame` fixtures are created inline.

There are focused test files for coverage, feature, network, joint, depth-sparsity, scree, elbow, cross-talk, singleton-ratio, QC metrics, class conversion, and file loading.

## Branch Workflow

Per `CONTRIBUTING.md`, this project uses a **develop branch workflow**:

- `master` is protected and production-ready.
- `develop` is the integration branch; PRs target `develop`.
- Feature/fix/docs branches branch from `develop` using prefixes like `feature/`, `fix/`, `docs/`.

Commits follow conventional commits (`feat:`, `fix:`, `docs:`, `test:`, `refactor:`, `chore:`).

## Style Guidelines

- 4-space indentation, spaces around operators, opening brace on the same line.
- `snake_case` for functions and variables.
- 80-character line limit (100 for comments).
- All exported functions must have roxygen2 docs with `@param`, `@return`, and `@examples`.
- Use `\dontrun{}` for examples requiring external data or long runtime.
