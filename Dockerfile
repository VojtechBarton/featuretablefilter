# FeatureTableFilter Docker Image
# A comprehensive environment for microbiome feature table filtering

FROM rocker/r-ver:4.4.2

LABEL maintainer="Vojtech Barton <vojtech.barton@hotmail.com>"
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
RUN Rscript -e "install.packages('remotes', repos = 'https://cloud.r-project.org')"

RUN Rscript -e "remotes::install_cran(c('devtools', 'roxygen2', 'testthat', 'knitr', 'rmarkdown', 'ggplot2', 'tidyr', 'dplyr', 'purrr', 'scales', 'patchwork', 'pheatmap', 'vegan', 'igraph', 'zoo', 'covr'), upgrade = 'never', Ncpus = 4)"

# Install Bioconductor packages using BiocManager
# R 4.4.x corresponds to Bioconductor 3.19; specify version to avoid conflicts
RUN Rscript -e "install.packages('BiocManager', repos = 'https://cran.r-project.org'); \
    BiocManager::install(version = '3.19', ask = FALSE, update = TRUE); \
    BiocManager::install(c('S4Vectors', 'SummarizedExperiment', 'SingleCellExperiment', 'TreeSummarizedExperiment', 'phyloseq'), ask = FALSE, update = FALSE, Ncpus = 4)"

# Copy package source to container
COPY . /workspace/featuretablefilter/

# Generate documentation and NAMESPACE with roxygen2 before installation
RUN Rscript -e "roxygen2::roxygenise('/workspace/featuretablefilter')"

# Install the featuretablefilter package from local source using Rscript
# This properly checks library paths for dependencies
RUN Rscript -e "install.packages('/workspace/featuretablefilter', repos = NULL, type = 'source', dependencies = FALSE)"
RUN echo "featuretablefilter installed successfully!"

# Verify installation
RUN Rscript -e "library(featuretablefilter); cat('\n=== Installed Functions ===\n'); funcs <- ls(asNamespace('featuretablefilter')); exported <- funcs[!startsWith(funcs, '.')]; cat(paste(sort(exported), collapse = '\n'), '\n'); cat('\n=== Package Version ===\n'); packageVersion('featuretablefilter')"

# Set working directory for user code
WORKDIR /workspace

# Default command
CMD ["R", "--vanilla"]
