#' Filter features using relative abundance cutoff applied as absolute threshold
#'
#' Calculates an absolute abundance threshold based on a relative percentage
#' of the least covered sample, then applies this threshold to filter features.
#' This ensures consistent filtering across all samples.
#'
#' @param table A feature table (data.frame or matrix) with features as rows and samples as columns.
#'              First column should be feature IDs, remaining columns are sample counts.
#' @param min_coverage Minimum total reads required per sample. Samples below this are removed first.
#' @param relative_threshold Relative abundance threshold as a proportion (e.g., 0.01 for 1\%).
#'                           This is converted to an absolute threshold based on the
#'                           minimum coverage sample.
#' @param remove_features Logical. If TRUE, remove features that fall below the threshold
#'                        in all samples. Default is TRUE.
#'
#' @return A list containing:
#' \describe{
#'   \item{table}{The filtered feature table}
#'   \item{absolute_threshold}{The calculated absolute threshold (relative_threshold * min_sample_coverage)}
#'   \item{min_sample_coverage}{Coverage of the least covered sample after initial filtering}
#'   \item{n_samples_removed}{Number of samples removed due to low coverage}
#'   \item{n_features_removed}{Number of features removed (if remove_features = TRUE)}
#' }
#'
#' @export
#'
#' @examples
#' data(example_feature_table)
#' result <- filter_by_relative_cutoff(
#'   example_feature_table,
#'   min_coverage = 1000,
#'   relative_threshold = 0.01
#' )
#' nrow(result$table)
filter_by_relative_cutoff <- function(table, min_coverage, relative_threshold,
                                       remove_features = TRUE) {
  # Step 1: Filter samples by minimum coverage
  sample_sums <- colSums(table[, -1, drop = FALSE])
  keep_samples <- sample_sums >= min_coverage

  n_samples_removed <- sum(!keep_samples)
  table_filtered <- table[, c(TRUE, keep_samples)]

  # Step 2: Find the minimum coverage among remaining samples
  sample_sums_filtered <- colSums(table_filtered[, -1, drop = FALSE])
  min_sample_coverage <- min(sample_sums_filtered)

  # Step 3: Calculate absolute threshold based on relative threshold
  absolute_threshold <- relative_threshold * min_sample_coverage

  # Step 4: Filter features below absolute threshold
  abundances <- table_filtered[, -1, drop = FALSE]
  feature_totals <- rowSums(abundances)

  if (remove_features) {
    # Remove features with total abundance below threshold
    keep_features <- feature_totals >= absolute_threshold
    table_final <- table_filtered[keep_features, , drop = FALSE]
    n_features_removed <- sum(!keep_features)
  } else {
    # Just zero out low-abundance features
    abundances[abundances < absolute_threshold] <- 0
    table_final <- cbind(table_filtered[, 1, drop = FALSE], abundances)
    n_features_removed <- 0
  }

  # Ensure column names are preserved after filtering
  if (is.matrix(table_final)) {
    table_final <- as.data.frame(table_final)
  }
  colnames(table_final) <- colnames(table_filtered)

  return(list(
    table = table_final,
    absolute_threshold = absolute_threshold,
    min_sample_coverage = min_sample_coverage,
    n_samples_removed = n_samples_removed,
    n_features_removed = n_features_removed
  ))
}
