\name{estimate_good_coverage}
\alias{estimate_good_coverage}
\title{Estimate Good's coverage for samples}
\description{
Calculates Good's coverage estimator, which measures the completeness of sampling
by estimating the probability that the next read will belong to a previously observed feature.
This is particularly useful for determining adequate sequencing depth in microbiome studies.

Good's coverage formula: C = 1 - (n_1 / n)
where n_1 is the number of singletons (features with exactly one read) and n is total reads.
}
\usage{
estimate_good_coverage(table, target_coverage = 0.95)
}
\arguments{
  \item{table}{A feature table (data.frame or matrix) with features as rows and samples as columns.
              First column should be feature IDs, remaining columns are sample counts.}
  \item{target_coverage}{Target coverage threshold (0-1). Samples below this coverage will be flagged.
                        Default is 0.95 (95\% coverage).}
}
\value{
A list containing:
  \item{sample_coverage}{Named vector of coverage estimates for each sample}
  \item{mean_coverage}{Mean coverage across all samples}
  \item{min_coverage}{Minimum coverage among all samples}
  \item{max_coverage}{Maximum coverage among all samples}
  \item{total_singletons}{Total number of singletons across all samples}
  \item{total_reads}{Total number of reads across all samples}
  \item{n_samples}{Total number of samples}
  \item{n_samples_below_target}{Number of samples below target coverage}
  \item{target_coverage}{The target coverage value used}
  \item{suggested_cutoff}{Suggested minimum read cutoff based on coverage analysis}
}
\examples{
# Estimate Good's coverage for all samples
# coverage_est <- estimate_good_coverage(table, target_coverage = 0.95)
# print(coverage_est)

# Filter samples to achieve at least 90\% coverage
# coverage_est <- estimate_good_coverage(table, target_coverage = 0.90)
}
