#' Estimate minimum read cutoff using Interquartile Range (IQR) / Tukey's Method
#'
#' Calculates a data-driven threshold for filtering low-coverage samples
#' based on the first quartile and IQR of sample coverages.
#'
#' @param table A feature table (data.frame or matrix) with features as rows and samples as columns.
#'              First column should be feature IDs, remaining columns are sample counts.
#' @param multiplier A numeric multiplier for the IQR. Samples with coverage below
#'                   (Q1 - multiplier * IQR) are considered outliers. Default is 1.5.
#' @param floor Minimum possible cutoff value. The returned cutoff will not be lower than this.
#'
#' @return A list containing:
#'   \item{cutoff}{The estimated minimum read threshold}
#'   \item{q1}{First quartile (25th percentile) of coverage}
#'   \item{q3}{Third quartile (75th percentile) of coverage}
#'   \item{iqr}{Interquartile range (Q3 - Q1)}
#'   \item{n_samples}{Total number of samples}
#'   \item{n_filtered}{Number of samples that would be filtered at this cutoff}
#'
#' @export
#'
#' @examples
#' # estimate_iqr_cutoff(my_table, multiplier = 1.5)
estimate_iqr_cutoff <- function(table, multiplier = 1.5, floor = 0) {
  # Calculate total reads per sample
  sample_sums <- colSums(table[, -1, drop = FALSE])

  # Calculate quartiles and IQR
  q1 <- quantile(sample_sums, 0.25)
  q3 <- quantile(sample_sums, 0.75)
  iqr <- q3 - q1

  # Calculate cutoff using Tukey's method (lower fence)
  cutoff <- q1 - multiplier * iqr

  # Apply floor
  cutoff <- max(cutoff, floor)

  # Count samples that would be filtered
  n_filtered <- sum(sample_sums < cutoff)

  return(list(
    cutoff = cutoff,
    q1 = q1,
    q3 = q3,
    iqr = iqr,
    n_samples = length(sample_sums),
    n_filtered = n_filtered
  ))
}
