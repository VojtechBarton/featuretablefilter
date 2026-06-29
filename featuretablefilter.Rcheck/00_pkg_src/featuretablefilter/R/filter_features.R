#' Filter features based on abundance threshold
#'
#' Sets low-abundance features to zero and removes features that are all zeros.
#' Works with both absolute counts and relative abundances.
#'
#' @param table A feature table (data.frame or matrix) with features as rows and samples as columns.
#'              First column should be feature IDs, remaining columns are sample counts/abundances.
#' @param threshold Minimum abundance threshold. Features with values below this will be set to zero.
#' @param mode Filtering mode: "absolute" for raw counts, "relative" for proportions/percentages.
#'             Default is "absolute".
#' @param remove_zeros Logical. If TRUE (default), remove features that become all zeros after filtering.
#' @param min_samples Optional. Minimum number of samples where a feature must exceed the threshold
#'                    to be kept. Default is 1 (feature must exceed threshold in at least one sample).
#'
#' @return A filtered feature table (data.frame) with low-abundance features zeroed and/or removed.
#'
#' @export
#'
#' @examples
#' # Filter features with absolute count < 5
#' # filtered <- filter_features_by_abundance(table, threshold = 5, mode = "absolute")
#'
#' # Filter features with relative abundance < 0.001 (0.1%)
#' # filtered <- filter_features_by_abundance(table, threshold = 0.001, mode = "relative")
#'
#' # Keep only features present in at least 3 samples
#' # filtered <- filter_features_by_abundance(table, threshold = 1, min_samples = 3)
filter_features_by_abundance <- function(table, threshold, mode = c("absolute", "relative"),
                                          remove_zeros = TRUE, min_samples = 1) {
  # Validate mode
  mode <- match.arg(mode)

  # Extract feature IDs and abundance columns
  feature_ids <- table[, 1, drop = FALSE]
  abundances <- table[, -1, drop = FALSE]

  # Convert to relative if requested
  if (mode == "relative") {
    sample_totals <- colSums(abundances)
    abundances_rel <- sweep(as.matrix(abundances), 2, sample_totals, FUN = "/")
    abundances_rel[is.na(abundances_rel)] <- 0
    # Preserve column names after sweep
    colnames(abundances_rel) <- colnames(abundances)
    abundances_rel <- as.data.frame(abundances_rel)
  } else {
    abundances_rel <- abundances
  }

  # Set values below threshold to zero
  filtered <- abundances_rel
  filtered[filtered < threshold] <- 0

  # Count samples where each feature exceeds threshold
  n_exceeding <- rowSums(filtered > 0)

  # Determine which features to keep
  if (remove_zeros) {
    # Keep features that exceed threshold in at least min_samples
    keep_features <- n_exceeding >= min_samples
  } else {
    # Zero out but don't remove - keep all features
    keep_features <- rep(TRUE, nrow(filtered))
  }

  # Build result
  result <- cbind(feature_ids[keep_features, , drop = FALSE],
                  filtered[keep_features, , drop = FALSE])

  # Ensure column names are preserved
  if (is.matrix(result)) {
    result <- as.data.frame(result)
  }
  colnames(result) <- colnames(table)

  # Return attributes about what was filtered
  attr(result, "n_filtered_out") <- sum(!keep_features)
  attr(result, "threshold") <- threshold
  attr(result, "mode") <- mode

  return(result)
}
