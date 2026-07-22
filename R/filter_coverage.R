#' Filter samples by minimum coverage
#'
#' Removes samples with total read count below the threshold.
#'
#' @param table A feature table (data.frame or matrix) with features as rows and samples as columns.
#'              First column should be feature IDs, remaining columns are sample counts.
#' @param min_reads Minimum number of reads required to keep a sample.
#'
#' @return A filtered feature table (data.frame) with low-coverage samples removed.
#'
#' @export
#'
#' @examples
#' data(example_feature_table)
#' result <- filter_by_coverage(example_feature_table, min_reads = 1000)
#' ncol(result)
filter_by_coverage <- function(table, min_reads) {
  # Calculate total reads per sample (sum of all rows for each sample column)
  sample_sums <- colSums(table[, -1, drop = FALSE])

  # Find samples that meet the threshold
  keep_samples <- sample_sums >= min_reads

  # Keep first column (feature IDs) and samples meeting threshold
  result <- table[, c(TRUE, keep_samples), drop = FALSE]

  # Ensure column names are preserved
  if (is.matrix(result)) {
    result <- as.data.frame(result)
  }

  # Explicitly preserve original column names
  orig_colnames <- colnames(table)
  new_colnames <- orig_colnames[c(TRUE, keep_samples)]
  colnames(result) <- new_colnames

  return(result)
}
