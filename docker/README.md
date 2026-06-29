# Docker Setup for featuretablefilter

This directory contains Docker configuration files for running the featuretablefilter package in a containerized environment.

## Quick Start

### Build the Docker image

```bash
docker build -t featuretablefilter .
```

Or using docker-compose:

```bash
docker-compose build
```

### Run an interactive R session

```bash
# Using docker
docker run -it --rm -v $(pwd):/workspace featuretablefilter R

# Using docker-compose
docker-compose run r-session
```

### Run an R script

```bash
# Place your script in the current directory, then:
docker run -it --rm -v $(pwd):/workspace featuretablefilter Rscript your_script.R

# Or with docker-compose
docker-compose run run-script
```

### Run tests

```bash
docker-compose run test
```

## Available Services

| Service | Description |
|---------|-------------|
| `featuretablefilter` | General purpose service for mounting and working with files |
| `r-session` | Interactive R session with the package loaded |
| `run-script` | Run a specific R script |
| `test` | Run the package test suite |

## Included Packages

### Core Dependencies
- R 4.3.3 (from rocker/r-ver)
- ggplot2, tidyr, dplyr, purrr
- vegan, igraph, zoo
- pheatmap, patchwork, scales

### Bioconductor
- phyloseq
- SummarizedExperiment
- SingleCellExperiment
- TreeSummarizedExperiment

### Development Tools
- devtools, roxygen2
- testthat
- knitr, rmarkdown
- covr (for coverage testing)

## Example Usage

### Analyze a feature table

```bash
# Start an interactive session
docker-compose run r-session

# Inside R
library(featuretablefilter)
table <- load_feature_table("/workspace/my-data/feature-table.tsv")
filtered <- filter_by_coverage(table, min_reads = 1000)
```

### Run a batch analysis

```bash
# Create your analysis script
cat > analysis.R << 'EOF'
library(featuretablefilter)

# Load data
table <- load_feature_table("data/feature-table.tsv")

# Filter
filtered <- filter_by_network_connectivity(table)

# Save results
write.table(filtered, "results/filtered-table.tsv", sep="\t", row.names=FALSE)
EOF

# Run the script
docker-compose run run-script analysis.R
```

## Volume Mounting

The Docker setup mounts the current directory to `/workspace` inside the container. This allows you to:

1. Access your local data files from within the container
2. Save output files that appear in your local directory

Example directory structure:

```
your-project/
├── Dockerfile
├── docker-compose.yml
├── data/
│   └── feature-table.tsv
├── results/
│   └── (output files will appear here)
└── scripts/
    └── analysis.R
```

## Building for Production

For a smaller production image, you can use a multi-stage build. Add this to your Dockerfile:

```dockerfile
# Production-ready slim image
FROM featuretablefilter:latest as builder

# Install package only (no dev tools)
FROM rocker/r-base:4.3.3
COPY --from=builder /usr/local/lib/R/site-library/featuretablefilter /usr/local/lib/R/site-library/featuretablefilter
```

## Troubleshooting

### Permission issues with output files

If you encounter permission issues when writing files from the container:

```bash
# Run with user ID matching your host user
docker run -u $(id -u):$(id -g) -v $(pwd):/workspace featuretablefilter ...
```

### Slow build times

Use build cache effectively:

```bash
# Use BuildKit for faster builds
DOCKER_BUILDKIT=1 docker build -t featuretablefilter .
```

## License

Same as the featuretablefilter package (GPL-3)
