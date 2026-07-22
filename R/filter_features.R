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
#' data(example_feature_table)
#' result <- filter_features_by_abundance(example_feature_table, threshold = 10)
#' nrow(result)
filter_features_by_abundance <- function(table, threshold, mode = c("absolute", "relative"),
                                          remove_zeros = TRUE, min_samples = 1) {
  # Validate mode
  mode <- match.arg(mode)

  # Extract feature IDs and abundance columns (keep original for output)
  feature_ids <- table[, 1, drop = FALSE]
  abundances_orig <- table[, -1, drop = FALSE]
  abundances <- as.matrix(abundances_orig)

  # Convert to relative if requested (for decision making only)
  if (mode == "relative") {
    sample_totals <- colSums(abundances)
    abundances_rel <- sweep(abundances, 2, sample_totals, FUN = "/")
    abundances_rel[is.na(abundances_rel)] <- 0
  } else {
    abundances_rel <- abundances
  }

  # Set values below threshold to zero (using relative for decision)
  mask <- abundances_rel < threshold
  n_exceeding <- rowSums(!mask)

  # Determine which features to keep
  if (remove_zeros) {
    # Keep features that exceed threshold in at least min_samples
    keep_features <- n_exceeding >= min_samples
  } else {
    # Zero out but don't remove - keep all features
    keep_features <- rep(TRUE, nrow(abundances))
  }

  # Apply mask to ORIGINAL abundances (not relative!)
  filtered <- abundances
  filtered[mask] <- 0

  # Build result using original count values
  result <- cbind(feature_ids[keep_features, , drop = FALSE],
                  filtered[keep_features, , drop = FALSE])

  # Ensure data.frame and column names are preserved
  result <- as.data.frame(result)
  colnames(result) <- colnames(table)

  # Return attributes about what was filtered
  attr(result, "n_filtered_out") <- sum(!keep_features)
  attr(result, "threshold") <- threshold
  attr(result, "mode") <- mode

  return(result)
}
