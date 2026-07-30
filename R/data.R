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
