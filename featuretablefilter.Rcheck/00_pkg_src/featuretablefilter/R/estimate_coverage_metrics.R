#' Estimate Good's coverage for samples
#'
#' Calculates Good's coverage estimator, which measures the completeness of sampling
#' by estimating the probability that the next read will belong to a previously observed feature.
#' This is particularly useful for determining adequate sequencing depth in microbiome studies.
#'
#' Good's coverage formula: C = 1 - (n_1 / n)
#' where n_1 is the number of singletons (features with exactly one read) and n is total reads.
#'
#' @param table A feature table (data.frame or matrix) with features as rows and samples as columns.
#'              First column should be feature IDs, remaining columns are sample counts.
#' @param target_coverage Target coverage threshold (0-1). Samples below this coverage will be flagged.
#'                        Default is 0.95 (95% coverage).
#'
#' @return A list containing:
#'   \item{sample_coverage}{Named vector of coverage estimates for each sample}
#'   \item{mean_coverage}{Mean coverage across all samples}
#'   \item{min_coverage}{Minimum coverage among all samples}
#'   \item{max_coverage}{Maximum coverage among all samples}
#'   \item{total_singletons}{Total number of singletons across all samples}
#'   \item{total_reads}{Total number of reads across all samples}
#'   \item{n_samples}{Total number of samples}
#'   \item{n_samples_below_target}{Number of samples below target coverage}
#'   \item{target_coverage}{The target coverage value used}
#'   \item{suggested_cutoff}{Suggested minimum read cutoff based on coverage analysis}
#'
#' @export
#'
#' @examples
#' # Estimate Good's coverage for all samples
#' # coverage_est <- estimate_good_coverage(table, target_coverage = 0.95)
#' # print(coverage_est)
#'
#' # Filter samples to achieve at least 90% coverage
#' # coverage_est <- estimate_good_coverage(table, target_coverage = 0.90)
estimate_good_coverage <- function(table, target_coverage = 0.95) {
  # Extract abundance matrix (exclude feature ID column)
  abundances <- as.matrix(table[, -1, drop = FALSE])

  # Calculate total reads per sample
  sample_totals <- colSums(abundances)

  # Count singletons per sample (features with exactly 1 read)
  singletons_per_sample <- apply(abundances, 2, function(x) sum(x == 1))

  # Calculate Good's coverage for each sample
  # C = 1 - (singletons / total_reads)
  sample_coverage <- 1 - (singletons_per_sample / sample_totals)

  # Handle edge cases (samples with 0 reads)
  sample_coverage[sample_totals == 0] <- NA

  # Count samples below target coverage
  n_below_target <- sum(sample_coverage < target_coverage, na.rm = TRUE)

  # Determine suggested cutoff: find minimum reads needed to achieve target coverage
  # This is heuristic: suggest filtering out the lowest coverage samples
  suggested_cutoff <- NA
  if (n_below_target > 0) {
    # Sort samples by coverage and find the cutoff point
    sorted_coverage <- sort(sample_coverage, decreasing = TRUE, na.last = TRUE)
    # Suggest keeping samples up to the point where coverage stabilizes
    cumulative_reads <- cumsum(rev(sort(sample_totals)))
    total_reads_all <- sum(sample_totals)

    # Find point where we've captured most reads while maintaining good coverage
    high_cov_samples <- names(sorted_coverage)[!is.na(sorted_coverage) & sorted_coverage >= target_coverage]
    if (length(high_cov_samples) > 0) {
      suggested_cutoff <- min(sample_totals[names(sample_totals) %in% high_cov_samples])
    }
  }

  return(list(
    sample_coverage = sample_coverage,
    mean_coverage = mean(sample_coverage, na.rm = TRUE),
    min_coverage = min(sample_coverage, na.rm = TRUE),
    max_coverage = max(sample_coverage, na.rm = TRUE),
    total_singletons = sum(singletons_per_sample, na.rm = TRUE),
    total_reads = sum(sample_totals),
    n_samples = ncol(abundances),
    n_samples_below_target = n_below_target,
    target_coverage = target_coverage,
    suggested_cutoff = suggested_cutoff
  ))
}

#' Estimate Chao's sample coverage for samples
#'
#' Calculates Chao's sample coverage estimator, a more conservative measure of
#' sampling completeness that accounts for unseen species/features. This estimator
#' is based on the frequency of rare features (singletons and doubletons).
#'
#' Chao's coverage formula (simplified):
#' C_hat = 1 - (f1/S) + (f1/n) * ((n-1)/n) * ((S-1)/S)
#' where f1 = singletons, S = number of features, n = total reads
#'
#' @param table A feature table (data.frame or matrix) with features as rows and samples as columns.
#'              First column should be feature IDs, remaining columns are sample counts.
#' @param target_coverage Target coverage threshold (0-1). Samples below this coverage will be flagged.
#'                        Default is 0.90 (90% coverage).
#'
#' @return A list containing:
#'   \item{sample_coverage}{Named vector of Chao coverage estimates for each sample}
#'   \item{mean_coverage}{Mean coverage across all samples}
#'   \item{min_coverage}{Minimum coverage among all samples}
#'   \item{max_coverage}{Maximum coverage among all samples}
#'   \item{total_singletons}{Total number of singletons across all samples}
#'   \item{total_doubletons}{Total number of doubletons across all samples}
#'   \item{total_features}{Total number of features across all samples}
#'   \item{total_reads}{Total number of reads across all samples}
#'   \item{n_samples}{Total number of samples}
#'   \item{n_samples_below_target}{Number of samples below target coverage}
#'   \item{target_coverage}{The target coverage value used}
#'   \item{suggested_cutoff}{Suggested minimum read cutoff based on coverage analysis}
#'
#' @export
#'
#' @examples
#' # Estimate Chao's coverage for all samples
#' # coverage_est <- estimate_chao_coverage(table, target_coverage = 0.90)
#' # print(coverage_est)
#'
#' # Compare Good's vs Chao's coverage
#' # good_cov <- estimate_good_coverage(table, target_coverage = 0.95)
#' # chao_cov <- estimate_chao_coverage(table, target_coverage = 0.90)
estimate_chao_coverage <- function(table, target_coverage = 0.90) {
  # Extract abundance matrix (exclude feature ID column)
  abundances <- as.matrix(table[, -1, drop = FALSE])

  # Calculate totals per sample
  sample_totals <- colSums(abundances)
  n_features_per_sample <- colSums(abundances > 0)

  # Count singletons and doubletons per sample
  singletons_per_sample <- apply(abundances, 2, function(x) sum(x == 1))
  doubletons_per_sample <- apply(abundances, 2, function(x) sum(x == 2))

  # Calculate Chao's coverage for each sample
  # C_hat = 1 - (f1/S) + (f1/n) * ((n-1)/n) * ((S-1)/S)
  sample_coverage <- sapply(seq_along(sample_totals), function(i) {
    n <- sample_totals[i]
    S <- n_features_per_sample[i]
    f1 <- singletons_per_sample[i]

    if (n == 0 || S == 0) return(NA)

    # Chao's estimator
    coverage <- 1 - (f1 / S) + (f1 / n) * ((n - 1) / n) * ((S - 1) / S)

    # Ensure coverage is between 0 and 1
    coverage <- max(0, min(1, coverage))
    coverage
  })

  names(sample_coverage) <- colnames(abundances)

  # Count samples below target coverage
  n_below_target <- sum(sample_coverage < target_coverage, na.rm = TRUE)

  # Determine suggested cutoff
  suggested_cutoff <- NA
  if (n_below_target > 0) {
    sorted_coverage <- sort(sample_coverage, decreasing = TRUE, na.last = TRUE)
    high_cov_samples <- names(sorted_coverage)[!is.na(sorted_coverage) & sorted_coverage >= target_coverage]
    if (length(high_cov_samples) > 0) {
      suggested_cutoff <- min(sample_totals[names(sample_totals) %in% high_cov_samples])
    }
  }

  return(list(
    sample_coverage = sample_coverage,
    mean_coverage = mean(sample_coverage, na.rm = TRUE),
    min_coverage = min(sample_coverage, na.rm = TRUE),
    max_coverage = max(sample_coverage, na.rm = TRUE),
    total_singletons = sum(singletons_per_sample, na.rm = TRUE),
    total_doubletons = sum(doubletons_per_sample, na.rm = TRUE),
    total_features = sum(n_features_per_sample),
    total_reads = sum(sample_totals),
    n_samples = ncol(abundances),
    n_samples_below_target = n_below_target,
    target_coverage = target_coverage,
    suggested_cutoff = suggested_cutoff
  ))
}
