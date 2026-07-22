#' Example feature table
#'
#' A synthetic microbiome feature table containing 50 features (ASVs) across
#' 20 samples. Provided for testing and documentation.
#'
#' @format A data frame with 50 rows (features) and 21 columns. The first column
#'   (#OTU ID) gives the feature identifier (ASV_001 to ASV_050) and the remaining
#'   20 columns (Sample_001 to Sample_020) contain read counts for each sample.
#'
#' @source Synthetic data generated for package testing and demonstration.
#'
#' @examples
#' data(example_feature_table)
#' str(example_feature_table)
#' head(example_feature_table)
#'
#' @docType data
#' @keywords datasets
#' @name example_feature_table
NULL

#' Example phyloseq object
#'
#' A synthetic phyloseq object containing 100 features across 30 samples.
#' Provided for testing phyloseq-compatible functions.
#'
#' @format A phyloseq object with otu_table, tax_table, sam_data, and phy_tree.
#'
#' @source Synthetic data generated for package testing and demonstration.
#'   Located in `inst/extdata/` rather than `data/` because it requires the
#'   \code{phyloseq} package to load.
#'
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   example_phyloseq <- readRDS(system.file("extdata", "example_phyloseq_object.rds",
#'     package = "featuretablefilter"))
#'   print(example_phyloseq)
#' }
#'
#' @seealso \code{link[phyloseq]{phyloseq}}
#'
#' @name example_phyloseq
#' @keywords internal
NULL

#' Example TreeSummarizedExperiment object
#'
#' A synthetic TreeSummarizedExperiment object containing 100 features across
#' 30 samples. Provided for testing TreeSummarizedExperiment-compatible
#' functions.
#'
#' @format A TreeSummarizedExperiment object with rowData, colData, and
#'   rowTree.
#'
#' @source Synthetic data generated for package testing and demonstration.
#'   Located in `inst/extdata/` rather than `data/` because it requires the
#'   \code{TreeSummarizedExperiment} package to load.
#'
#' @examples
#' if (requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
#'   example_treesummarizedexperiment <- readRDS(system.file("extdata",
#'     "example_treesummarizedexperiment_object.rds", package = "featuretablefilter"))
#'   print(example_treesummarizedexperiment)
#' }
#'
#' @seealso \code{link[TreeSummarizedExperiment]{TreeSummarizedExperiment}}
#'
#' @name example_treesummarizedexperiment
#' @keywords internal
NULL
