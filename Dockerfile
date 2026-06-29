# FeatureTableFilter Docker Image
# A comprehensive environment for microbiome feature table filtering

FROM rocker/r-ver:4.3.3

LABEL maintainer="Vojtech Barton <vojtech.barton@gmail.com>"
LABEL description="Docker image with featuretablefilter R package and all dependencies"

# Set working directory
WORKDIR /workspace

# Update system packages and install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgdal-dev \
    libproj-dev \
    libgeos-dev \
    libgit2-dev \
    libv8-dev \
    default-mysql-client \
    default-libmysqlclient-dev \
    build-essential \
    pkg-config \
    cmake \
    git \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

# Install R packages from CRAN
RUN Rscript -e 'install.packages("remotes", repos = "https://cloud.r-project.org")'

RUN Rscript -e 'remotes::install_cran(c(\
    "devtools",\
    "roxygen2",\
    "testthat",\
    "knitr",\
    "rmarkdown",\
    "ggplot2",\
    "tidyr",\
    "dplyr",\
    "purrr",\
    "scales",\
    "patchwork",\
    "pheatmap",\
    "vegan",\
    "igraph",\
    "zoo",\
    "covr"\
  ), \
  upgrade = "never", \
  Ncpus = 4)'

# Install Bioconductor packages
RUN Rscript -e '\
if (!requireNamespace("BiocManager", quietly = TRUE)) \
    install.packages("BiocManager", repos = "https://bioconductor.org/packages/3.18/bioc", upgrade = FALSE)\
BiocManager::install(version = "3.18", ask = FALSE, force = TRUE)\
remotes::install_bioc(c(\
    "SummarizedExperiment",\
    "SingleCellExperiment",\
    "TreeSummarizedExperiment",\
    "phyloseq"\
  ), \
  upgrade = "never", \
  Ncpus = 4)'

# Copy package source to container
COPY . /workspace/featuretablefilter/

# Install the featuretablefilter package from local source
RUN Rscript -e '\
  setwd("/workspace/featuretablefilter")\
  devtools::install(".", dependencies = FALSE, upgrade = "never")\
  cat("featuretablefilter installed successfully!\n")'

# Verify installation
RUN Rscript -e '\
  library(featuretablefilter)\
  cat("\n=== Installed Functions ===\n")\
  funcs <- ls(asNamespace("featuretablefilter"))\
  exported <- funcs[grepl("^[^\\.]", funcs)]\
  cat(paste(sort(exported), collapse = "\n"), "\n")\
  cat("\n=== Package Version ===\n")\
  packageVersion("featuretablefilter")\
'

# Set working directory for user code
WORKDIR /workspace

# Default command
CMD ["R", "--vanilla"]
