#' Estimate minimum read cutoff using Median Absolute Deviation (MAD)
#'
#' Calculates a data-driven threshold for filtering low-coverage samples
#' based on the median and MAD of sample coverages.
#'
#' @param table A feature table (data.frame or matrix) with features as rows and samples as columns.
#'              First column should be feature IDs, remaining columns are sample counts.
#' @param multiplier A numeric multiplier for the MAD. Samples with coverage below
#'                   (median - multiplier * MAD) are considered outliers. Default is 3.
#' @param floor Minimum possible cutoff value. The returned cutoff will not be lower than this.
#'
#' @return A list containing:
#'   \item{cutoff}{The estimated minimum read threshold}
#'   \item{median}{Median coverage across samples}
#'   \item{mad}{Median Absolute Deviation of coverage}
#'   \item{n_samples}{Total number of samples}
#'   \item{n_filtered}{Number of samples that would be filtered at this cutoff}
#'
#' @export
#'
#' @examples
#' # estimate_mad_cutoff(my_table, multiplier = 3)
estimate_mad_cutoff <- function(table, multiplier = 3, floor = 0) {
  # Calculate total reads per sample
  sample_sums <- colSums(table[, -1, drop = FALSE])

  # Calculate median and MAD
  median_cov <- median(sample_sums)
  mad_cov <- mad(sample_sums)

  # Calculate cutoff (lower bound for "normal" samples)
  cutoff <- median_cov - multiplier * mad_cov

  # Apply floor
  cutoff <- max(cutoff, floor)

  # Count samples that would be filtered
  n_filtered <- sum(sample_sums < cutoff)

  return(list(
    cutoff = cutoff,
    median = median_cov,
    mad = mad_cov,
    n_samples = length(sample_sums),
    n_filtered = n_filtered
  ))
}
