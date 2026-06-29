#' Filter samples by singleton/doubleton ratio
#'
#' Removes samples with an unnaturally high ratio of singleton (count=1) or
#' doubleton (count=2) features relative to their total sequencing depth.
#' High ratios often indicate poor sequencing runs or severe PCR artifacts.
#'
#' @param table A feature table (data.frame or matrix) with features as rows and samples as columns.
#'              First column should be feature IDs, remaining columns are sample counts.
#' @param max_singleton_ratio Maximum allowed ratio of (singletons + doubletons) / total_reads.
#'                            Samples exceeding this ratio will be removed. Default is 0.1 (10%).
#' @param count_type Type of low-count features to consider: "singleton" (count=1 only),
#'                   "doubleton" (count=2 only), or "both" (counts 1 and 2). Default is "both".
#'
#' @return A filtered feature table (data.frame) with suspicious samples removed.
#'         The result includes attributes:
#'         \itemize{
#'           \item \code{n_filtered_out}: Number of samples removed
#'           \item \code{n_retained}: Number of samples retained
#'           \item \code{threshold}: The max_singleton_ratio used
#'           \item \code{ratio_vector}: Named vector of singleton ratios for all original samples
#'         }
#'
#' @export
#'
#' @examples
#' # Remove samples where >10% of reads are in singletons/doubletons
#' # filtered_table <- filter_by_singleton_ratio(my_table, max_singleton_ratio = 0.1)
#'
#' # Only consider true singletons (count = 1)
#' # filtered_table <- filter_by_singleton_ratio(my_table, max_singleton_ratio = 0.05, count_type = "singleton")
filter_by_singleton_ratio <- function(table, max_singleton_ratio = 0.1, count_type = "both") {
  # Validate inputs
  if (!is.data.frame(table) && !is.matrix(table)) {
    stop("table must be a data.frame or matrix")
  }

  if (ncol(table) < 2) {
    stop("table must have at least one sample column (beyond feature ID column)")
  }

  if (!is.numeric(max_singleton_ratio) || max_singleton_ratio <= 0 || max_singleton_ratio > 1) {
    stop("max_singleton_ratio must be a numeric value between 0 and 1")
  }

  if (!count_type %in% c("singleton", "doubleton", "both")) {
    stop("count_type must be 'singleton', 'doubleton', or 'both'")
  }

  # Get sample columns (exclude first column which is feature IDs)
  sample_data <- table[, -1, drop = FALSE]

  # Calculate total reads per sample
  sample_totals <- colSums(sample_data)

  # Count singletons (features with count = 1) per sample
  singleton_counts <- colSums(sample_data == 1)

  # Count doubletons (features with count = 2) per sample
  doubleton_counts <- colSums(sample_data == 2)

  # Determine which counts to use based on count_type
  if (count_type == "singleton") {
    low_count_sum <- singleton_counts
  } else if (count_type == "doubleton") {
    low_count_sum <- doubleton_counts
  } else {  # "both"
    low_count_sum <- singleton_counts + doubleton_counts
  }

  # Calculate ratio of low-count reads to total reads per sample
  # Handle edge case of zero total reads
  ratio_vector <- ifelse(
    sample_totals > 0,
    low_count_sum / sample_totals,
    NA
  )

  # Name the ratio vector
  names(ratio_vector) <- colnames(sample_data)

  # Identify samples to keep (ratio below threshold or NA due to zero reads)
  keep_samples <- is.na(ratio_vector) | (ratio_vector <= max_singleton_ratio)

  # Create filtered table
  result <- table[, c(TRUE, keep_samples), drop = FALSE]

  # Ensure data.frame structure is preserved
  if (is.matrix(result)) {
    result <- as.data.frame(result)
  }

  # Preserve column names
  orig_colnames <- colnames(table)
  new_colnames <- orig_colnames[c(TRUE, keep_samples)]
  colnames(result) <- new_colnames

  # Attach informative attributes
  attr(result, "n_filtered_out") <- sum(!keep_samples)
  attr(result, "n_retained") <- sum(keep_samples)
  attr(result, "threshold") <- max_singleton_ratio
  attr(result, "ratio_vector") <- ratio_vector

  return(result)
}
