# Shiny Dashboard Plan

## Overview

Interactive web application for exploring filtering options and visualizing results in real-time.

## Features

### 1. Data Upload
- File upload (TSV, CSV)
- Support for file path input (server deployment)
- Auto-detect format and display summary statistics

### 2. Interactive Parameter Selection

#### Coverage Filtering
- Method selector: none, absolute, mad, iqr, good, chao
- Dynamic parameter inputs based on method:
  - `absolute`: min_reads slider/input
  - `mad`/`iqr`: multiplier slider, floor input
  - `good`/`chao`: target coverage slider (0-1)

#### Singleton Ratio Filtering
- Method selector: none, absolute
- max_singleton_ratio slider
- count_type radio buttons: singleton, doubleton, both

#### Cross-Talk Filtering
- Method selector: none, zero, remove_feature, flag
- max_rel_threshold slider
- min_abs_cutoff input (optional)

#### Abundance Filtering
- Method selector: none, absolute, relative, relative_cutoff, joint
- Dynamic parameters based on method
- For joint: AND/OR logic selector

### 3. Real-time Preview
- Updated table dimensions after each parameter change
- Retention rates (features, samples, reads)
- Sparsity changes
- Processing status indicator

### 4. Visualizations
- Coverage histogram with threshold line
- Feature abundance distribution
- Sparsity comparison (before/after)
- Retention rates barplot
- Scree plot (if enabled)
- All plots update dynamically

### 5. Results Export
- Download filtered table (TSV)
- Download QC report (text)
- Download all plots (ZIP)
- Export R code for reproducibility

## File Structure

```
inst/shiny-dashboard/
├── app.R                 # Main Shiny application
├── server.R              # Server logic
├── ui.R                  # User interface
└── www/                  # Static files (CSS, JS)
    └── style.css
```

## Dependencies

```r
# In DESCRIPTION
Suggests:
    shiny,
    DT,
    downloadhelper
```

## Usage

```r
# From package
library(featuretablefilter)
runDashboard()

# Direct
shiny::runApp(system.file("shiny-dashboard", package = "featuretablefilter"))
```

## Development Tasks

- [ ] Create directory structure
- [ ] Design UI layout (sidebar + main panel)
- [ ] Implement data upload handler
- [ ] Create reactive parameter inputs
- [ ] Build server-side filtering logic
- [ ] Implement real-time preview
- [ ] Add dynamic visualizations
- [ ] Create export functions
- [ ] Add loading indicators
- [ ] Test with various data sizes
- [ ] Add help tooltips
- [ ] Write documentation
