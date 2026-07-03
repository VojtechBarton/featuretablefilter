#' Filter cross-talk/index-switching artifacts
#'
#' Removes likely index-hopping or cross-contamination artifacts by comparing
#' each feature's abundance in each sample to its maximum abundance across all samples.
#' Reads appearing at very low levels relative to a feature's peak abundance are
#' likely leakage from highly abundant samples on the same flow cell.
#'
#' This filter is particularly useful for:
#' - Illumina patterned flow cells (ExAmp chemistry) where index hopping occurs
#' - Multiplexed runs with highly uneven library concentrations
#' - Detecting and removing "ghost" ASVs that appear sporadically at low counts
#'
#' @param table A feature table (data.frame or matrix) with features as rows and samples as columns.
#'              First column should be feature IDs, remaining columns are sample counts.
#' @param max_rel_threshold Maximum allowed relative abundance as a fraction of the feature's
#'                          maximum abundance across all samples. Values below this threshold
#'                          relative to the max are considered leakage and set to zero.
#'                          Default is 0.001 (0.1% of max).
#' @param min_abs_cutoff Optional. Minimum absolute count required to be retained regardless
#'                       of the relative threshold. Features with counts >= this value are kept.
#'                       Set to NULL (default) to disable. Useful for preserving any real
#'                       low-abundance presence (e.g., min_abs_cutoff = 2 keeps doubletons).
#' @param mode Control how filtering is applied: "zero" sets suspected leakage reads to zero (default),
#'             "remove_feature" removes entire feature if any leakage detected,
#'             "flag" adds attribute with leakage info but doesn't modify data.
#' @param return_details Logical. If TRUE, add attributes with detailed information about
#'                       which values were flagged as leakage. Default is FALSE.
#'
#' @return A filtered feature table (data.frame) with suspected cross-talk artifacts removed.
#'         The result includes attributes: n_leakage_zeros (number of individual cells set to zero),
#'         n_features_affected (number of features with at least one leakage zero), threshold
#'         (the max_rel_threshold used), min_abs_cutoff, and mode. If return_details = TRUE, also
#'         includes leakage_matrix (logical matrix indicating which cells were flagged) and
#'         feature_max (named vector of each feature's maximum abundance).
#'
#' @export
#'
#' @examples
#' # Remove reads < 0.1% of each feature's maximum abundance
#' # cleaned <- filter_cross_talk(my_table, max_rel_threshold = 0.001)
#'
#' # Stricter: 0.01% threshold with minimum absolute cutoff of 3 reads
#' # cleaned <- filter_cross_talk(my_table, max_rel_threshold = 0.0001, min_abs_cutoff = 3)
#'
#' # Get detailed leakage information
#' # result <- filter_cross_talk(my_table, return_details = TRUE)
#' # leakage_pattern <- attr(result, "leakage_matrix")
filter_cross_talk <- function(table, max_rel_threshold = 0.001, min_abs_cutoff = NULL,
                               mode = c("zero", "remove_feature", "flag"),
                               return_details = FALSE) {
  # Validate inputs
  if (!is.data.frame(table) && !is.matrix(table)) {
    stop("table must be a data.frame or matrix")
  }

  if (ncol(table) < 2) {
    stop("table must have at least one sample column (beyond feature ID column)")
  }

  if (!is.numeric(max_rel_threshold) || max_rel_threshold <= 0 || max_rel_threshold > 1) {
    stop("max_rel_threshold must be a numeric value between 0 and 1")
  }

  if (!is.null(min_abs_cutoff)) {
    if (!is.numeric(min_abs_cutoff) || min_abs_cutoff < 0) {
      stop("min_abs_cutoff must be a non-negative numeric value")
    }
  }

  mode <- match.arg(mode)

  # Extract feature IDs and abundance columns
  feature_ids <- table[, 1, drop = FALSE]
  col_names <- colnames(table)
  abundances <- as.matrix(table[, -1, drop = FALSE])

  # Calculate maximum abundance for each feature across all samples
  feature_max <- apply(abundances, 1, max)
  names(feature_max) <- feature_ids[, 1]

  # Handle edge case: features with zero max (all zeros)
  feature_max[is.na(feature_max) | feature_max == 0] <- 0

  # Calculate dynamic threshold for each feature (max * ratio)
  # This creates a matrix where each row has the same threshold value repeated
  dynamic_threshold <- outer(feature_max, rep(1, ncol(abundances))) * max_rel_threshold

  # Identify leakage: values > 0 AND < dynamic_threshold AND (if min_abs_cutoff set) < min_abs_cutoff
  leakage_mask <- abundances > 0 & abundances < dynamic_threshold

  if (!is.null(min_abs_cutoff)) {
    # Also require being below absolute cutoff
    leakage_mask <- leakage_mask & abundances < min_abs_cutoff
  }

  # Count leakage statistics
  n_leakage_zeros <- sum(leakage_mask)
  n_features_affected <- sum(rowSums(leakage_mask) > 0)

  # Apply filtering based on mode
  if (mode == "flag") {
    # Don't modify data, just return original
    result <- as.data.frame(cbind(feature_ids, abundances))
    colnames(result) <- col_names
    rownames(result) <- feature_ids[, 1]
  } else if (mode == "zero") {
    # Set leakage values to zero
    filtered_abundances <- abundances
    filtered_abundances[leakage_mask] <- 0
    result <- as.data.frame(cbind(feature_ids, filtered_abundances))
    colnames(result) <- col_names
    rownames(result) <- feature_ids[, 1]
  } else if (mode == "remove_feature") {
    # Remove any feature with at least one leakage reading
    keep_features <- rowSums(leakage_mask) == 0
    result <- as.data.frame(cbind(feature_ids[keep_features, , drop = FALSE],
                                   abundances[keep_features, , drop = FALSE]))
    colnames(result) <- col_names
    rownames(result) <- feature_ids[keep_features, 1]
    n_features_affected <- sum(!keep_features)
  }

  # Build attributes
  attr(result, "n_leakage_zeros") <- ifelse(mode == "remove_feature" || mode == "flag", NA_integer_, n_leakage_zeros)
  attr(result, "n_features_affected") <- n_features_affected
  attr(result, "threshold") <- max_rel_threshold
  attr(result, "min_abs_cutoff") <- min_abs_cutoff
  attr(result, "mode") <- mode

  # Add detailed info if requested
  if (return_details) {
    leakage_mask_df <- as.data.frame(leakage_mask)
    rownames(leakage_mask_df) <- feature_ids[, 1]
    attr(result, "leakage_matrix") <- as.matrix(leakage_mask_df)
    attr(result, "feature_max") <- feature_max
  }

  return(result)
}


#' Alias for filter_cross_talk with more intuitive name
#'
#' @rdname filter_cross_talk
#' @export
filter_index_hopping <- filter_cross_talk
