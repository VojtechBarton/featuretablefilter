#' Filter samples based on ecological coverage estimators
#'
#' Filters samples based on Good's or Chao's coverage estimates, which measure
#' the completeness of sampling rather than just absolute sequencing depth.
#' This approach removes samples that are likely undersampled ecologically.
#'
#' @param table A feature table (data.frame or matrix) with features as rows and samples as columns.
#'              First column should be feature IDs, remaining columns are sample counts.
#' @param method Coverage estimation method: "good" for Good's coverage, "chao" for Chao's coverage.
#' @param target_coverage Target coverage threshold (0-1). Samples below this will be filtered.
#'                        Default is 0.95 for Good's, 0.90 for Chao's.
#' @param min_reads Optional minimum absolute read count cutoff applied in addition to coverage.
#'                  Useful as a safety floor. Default is 0 (no additional cutoff).
#'
#' @return A list containing:
#'   \item{table}{Filtered feature table}
#'   \item{coverage_before}{Coverage estimates before filtering}
#'   \item{coverage_after}{Coverage estimates after filtering}
#'   \item{mean_coverage_before}{Mean coverage before filtering}
#'   \item{mean_coverage_after}{Mean coverage after filtering}
#'   \item{n_samples_before}{Number of samples before filtering}
#'   \item{n_samples_after}{Number of samples after filtering}
#'   \item{n_samples_filtered}{Number of samples removed}
#'   \item{method}{Method used ("good" or "chao")}
#'   \item{target_coverage}{Target coverage threshold used}
#'
#' @export
#'
#' @examples
#' # Filter using Good's coverage (keep samples with >= 95% coverage)
#' # result <- filter_by_coverage_estimator(table, method = "good", target_coverage = 0.95)
#'
#' # Filter using Chao's coverage (keep samples with >= 90% coverage)
#' # result <- filter_by_coverage_estimator(table, method = "chao", target_coverage = 0.90)
filter_by_coverage_estimator <- function(table, method = c("good", "chao"),
                                          target_coverage = NULL, min_reads = 0) {
  # Validate method
  method <- match.arg(method)

  # Set default target coverage if not specified
  if (is.null(target_coverage)) {
    target_coverage <- if (method == "good") 0.95 else 0.90
  }

  # Estimate coverage
  if (method == "good") {
    coverage_est <- estimate_good_coverage(table, target_coverage = target_coverage)
  } else {
    coverage_est <- estimate_chao_coverage(table, target_coverage = target_coverage)
  }

  # Store coverage before filtering
  coverage_before <- coverage_est$sample_coverage
  n_samples_before <- length(coverage_before)

  # Determine which samples to keep
  # Keep samples that meet BOTH the coverage threshold AND minimum reads requirement
  sample_totals <- colSums(as.matrix(table[, -1, drop = FALSE]))
  keep_samples <- (coverage_before >= target_coverage) & (sample_totals >= min_reads)

  # Also ensure we don't remove all samples - keep at least the best ones
  if (all(!keep_samples)) {
    # Keep the sample with highest coverage
    best_sample <- which.max(coverage_before)
    keep_samples <- rep(FALSE, length(keep_samples))
    keep_samples[best_sample] <- TRUE
  }

  # Filter the table
  filtered_table <- filter_by_coverage(table, min_reads = 0)  # Get all data first
  filtered_table <- filtered_table[, c(TRUE, keep_samples), drop = FALSE]

  # Preserve column names
  if (is.matrix(filtered_table)) {
    filtered_table <- as.data.frame(filtered_table)
  }
  orig_colnames <- colnames(table)
  new_colnames <- orig_colnames[c(TRUE, keep_samples)]
  colnames(filtered_table) <- new_colnames

  # Calculate coverage after filtering
  if (ncol(filtered_table) > 1) {
    if (method == "good") {
      coverage_after_est <- estimate_good_coverage(filtered_table, target_coverage = target_coverage)
    } else {
      coverage_after_est <- estimate_chao_coverage(filtered_table, target_coverage = target_coverage)
    }
    coverage_after <- coverage_after_est$sample_coverage
    mean_coverage_after <- coverage_after_est$mean_coverage
  } else {
    coverage_after <- numeric(0)
    mean_coverage_after <- NA
  }

  return(list(
    table = filtered_table,
    coverage_before = coverage_before,
    coverage_after = coverage_after,
    mean_coverage_before = mean(coverage_before, na.rm = TRUE),
    mean_coverage_after = mean_coverage_after,
    n_samples_before = n_samples_before,
    n_samples_after = ncol(filtered_table) - 1,
    n_samples_filtered = sum(!keep_samples),
    method = method,
    target_coverage = target_coverage
  ))
}
