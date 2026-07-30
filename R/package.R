#' @importFrom S4Vectors DataFrame
#' @importFrom zoo rollapply
#' @import ggplot2
#' @importFrom grDevices colorRampPalette dev.off png
#' @importFrom graphics boxplot hist
#' @importFrom stats cor cor.test dist lm mad median quantile sd setNames var
#' @importFrom utils head read.table tail write.table
#' @importFrom methods new
NULL

# Suppress R CMD check NOTE for ggplot2 aes() column references and
# other non-standard evaluation variables
utils::globalVariables(c(
  "..fill.color..", "abundance", "category", "collapse_rate", "count",
  "coverage", "degree", "depth", "derivative", "feature", "group",
  "label", "mean_sparsity", "metric", "mid", "outlier_type",
  "pct_features_retained", "pct_reads_retained", "pct_samples_retained",
  "percentage", "proportion", "richness", "sample_idx", "sample_name",
  "smoothed_richness", "sparsity", "table_type", "threshold",
  "total_reads", "type", "value", "x", "y"
))

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
#'   \item{\code{example_feature_table}}{Synthetic feature table (50 features x 20 samples),
#'     available via \code{data(example_feature_table)}.}
#'   \item{\code{example_phyloseq_object.rds}}{Synthetic phyloseq object in
#'     \code{inst/extdata/}, loaded via \code{readRDS(system.file("extdata",
#'     "example_phyloseq_object.rds", package = "featuretablefilter"))}.}
#'   \item{\code{example_treesummarizedexperiment_object.rds}}{Synthetic TSE object in
#'     \code{inst/extdata/}, loaded via \code{readRDS(system.file("extdata",
#'     "example_treesummarizedexperiment_object.rds", package = "featuretablefilter"))}.}
#' }
#'
#' @docType package
#' @name featuretablefilter
#' @keywords internal
"_PACKAGE"
