#' Filter features using joint abundance and prevalence criteria
#'
#' Combines abundance thresholds with prevalence requirements using AND/OR logic.
#' This allows flexible filtering where features can be retained based on being
#' either abundant, widespread, or both.
#'
#' Common use cases:
#' - AND logic: Keep only features that are BOTH abundant enough AND widely distributed
#' - OR logic: Keep features that are EITHER abundant enough OR widely distributed
#'   (useful for retaining rare but high-abundance features, or common but low-abundance features)
#'
#' @param table A feature table (data.frame or matrix) with features as rows and samples as columns.
#'              First column should be feature IDs, remaining columns are sample counts/abundances.
#' @param abundance_threshold Minimum abundance threshold. Features below this in a sample are considered "absent"
#'                            for the abundance criterion. Default is 0 (any non-zero value counts).
#' @param prevalence_threshold Proportion of samples (0-1) where feature must meet abundance threshold.
#'                             For example, 0.3 means the feature must exceed abundance_threshold in at least 30% of samples.
#' @param mode Abundance measurement mode: "absolute" for raw counts, "relative" for proportions.
#'             Default is "relative".
#' @param logic Logical operator: "AND" requires both criteria to be met, "OR" requires either criterion.
#'              Default is "OR".
#' @param remove_zeros Logical. If TRUE (default), remove features that don't meet the filtering criteria.
#'
#' @return A list containing:
#'   \item{table}{Filtered feature table (data.frame)}
#'   \item{n_features_before}{Number of features before filtering}
#'   \item{n_features_after}{Number of features after filtering}
#'   \item{n_features_removed}{Number of features removed}
#'   \item{n_by_abundance_only}{Features meeting only abundance criterion (not prevalence)}
#'   \item{n_by_prevalence_only}{Features meeting only prevalence criterion (not abundance)}
#'   \item{n_by_both}{Features meeting both criteria}
#'   \item{n_by_neither}{Features meeting neither criterion (these were removed)}
#'   \item{abundance_threshold}{The abundance threshold used}
#'   \item{prevalence_threshold}{The prevalence threshold used}
#'   \item{logic}{The logic operator used ("AND" or "OR")}
#'   \item{mode}{The abundance mode used ("absolute" or "relative")}
#'   \item{feature_details}{Data.frame with per-feature filtering details}
#'
#' @export
#'
#' @examples
#' # Keep features that are EITHER >0.1% relative abundance OR present in >=30% of samples
#' # result <- filter_features_joint(table,
#' #                                 abundance_threshold = 0.001,
#' #                                 prevalence_threshold = 0.3,
#' #                                 mode = "relative",
#' #                                 logic = "OR")
#'
#' # Keep features that are BOTH >5 reads AND present in >=50% of samples
#' # result <- filter_features_joint(table,
#' #                                 abundance_threshold = 5,
#' #                                 prevalence_threshold = 0.5,
#' #                                 mode = "absolute",
#' #                                 logic = "AND")
filter_features_joint <- function(table,
                                   abundance_threshold = 0,
                                   prevalence_threshold = 0.3,
                                   mode = c("relative", "absolute"),
                                   logic = c("OR", "AND"),
                                   remove_zeros = TRUE) {
  # Validate inputs
  mode <- match.arg(mode)
  logic <- match.arg(logic)

  if (abundance_threshold < 0) {
    stop("abundance_threshold must be non-negative")
  }
  if (prevalence_threshold < 0 || prevalence_threshold > 1) {
    stop("prevalence_threshold must be between 0 and 1")
  }

  # Extract feature IDs and abundance matrix (keep original for output)
  feature_ids <- table[, 1, drop = FALSE]
  abundances_orig <- as.matrix(table[, -1, drop = FALSE])
  abundances <- abundances_orig
  n_samples <- ncol(abundances)

  # Convert to relative if requested (for decision making only)
  if (mode == "relative") {
    sample_totals <- colSums(abundances)
    abundances_rel <- sweep(abundances, 2, sample_totals, FUN = "/")
    abundances_rel[is.na(abundances_rel)] <- 0
  } else {
    abundances_rel <- abundances
  }

  # Calculate which cells exceed abundance threshold (using relative if mode is relative)
  exceeds_abundance <- abundances_rel >= abundance_threshold

  # Calculate prevalence for each feature (proportion of samples exceeding threshold)
  feature_prevalence <- rowSums(exceeds_abundance) / n_samples

  # Determine which features meet each criterion
  meets_abundance <- rowSums(exceeds_abundance) > 0  # At least one sample exceeds threshold
  meets_prevalence <- feature_prevalence >= prevalence_threshold

  # Apply AND/OR logic
  if (logic == "AND") {
    keep_features <- meets_abundance & meets_prevalence
  } else {  # OR
    keep_features <- meets_abundance | meets_prevalence
  }

  # Calculate statistics
  n_by_abundance_only <- sum(meets_abundance & !meets_prevalence)
  n_by_prevalence_only <- sum(!meets_abundance & meets_prevalence)
  n_by_both <- sum(meets_abundance & meets_prevalence)
  n_by_neither <- sum(!meets_abundance & !meets_prevalence)

  # Apply mask to ORIGINAL abundances (not relative!)
  mask <- !exceeds_abundance & keep_features == FALSE
  filtered_abundances <- abundances_orig
  filtered_abundances[mask] <- 0

  # Build result table using original count values
  if (remove_zeros) {
    result_table <- cbind(
      feature_ids[keep_features, , drop = FALSE],
      filtered_abundances[keep_features, , drop = FALSE]
    )
  } else {
    result_table <- cbind(feature_ids, filtered_abundances)
  }

  # Ensure proper data.frame structure
  if (is.matrix(result_table)) {
    result_table <- as.data.frame(result_table)
  }
  colnames(result_table) <- colnames(table)

  # Create feature details dataframe (use relative values for display purposes)
  feature_details <- data.frame(
    feature_id = table[, 1],
    max_abundance = apply(abundances_rel, 1, max),
    mean_abundance = rowMeans(abundances_rel),
    prevalence = feature_prevalence,
    meets_abundance = meets_abundance,
    meets_prevalence = meets_prevalence,
    kept = keep_features,
    stringsAsFactors = FALSE
  )

  return(list(
    table = result_table,
    n_features_before = nrow(table),
    n_features_after = sum(keep_features),
    n_features_removed = sum(!keep_features),
    n_by_abundance_only = n_by_abundance_only,
    n_by_prevalence_only = n_by_prevalence_only,
    n_by_both = n_by_both,
    n_by_neither = n_by_neither,
    abundance_threshold = abundance_threshold,
    prevalence_threshold = prevalence_threshold,
    logic = logic,
    mode = mode,
    feature_details = feature_details
  ))
}
