#' @importFrom S4Vectors DataFrame
#' @importFrom zoo rollapply
#' @import ggplot2
NULL

# Null coalescing operator for default values
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' featuretablefilter: Feature Table Filtering for Microbiome Data
#'
#' \code{featuretablefilter} provides functions for filtering microbiome feature
#' tables based on coverage, abundance, and network connectivity criteria. The
#' package supports native handling of \code{phyloseq} and
#' \code{TreeSummarizedExperiment} objects, providing comprehensive quality control
#' metrics and visualizations for filtering decisions.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{run_filtering_pipeline}: Complete filtering workflow
#'   \item \code{filter_by_coverage}: Sample coverage filtering
#'   \item \code{filter_features_by_abundance}: Feature abundance filtering
#'   \item \code{compute_filtering_qc}: Calculate QC metrics
#'   \item \code{runDashboard}: Interactive Shiny dashboard
#' }
#'
#' @section Example Data:
#' \describe{
#'   \item{\code{example_feature_table}}{Synthetic feature table (50 features x 20 samples)}
#'   \item{\code{example_phyloseq}}{Synthetic phyloseq object (100 features x 30 samples)}
#'   \item{\code{example_treesummarizedexperiment}}{Synthetic TSE object (100 features x 30 samples)}
#' }
#'
#' @docType package
#' @name featuretablefilter
#' @keywords internal
"_PACKAGE"
