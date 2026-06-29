\name{filter_by_coverage_estimator}
\alias{filter_by_coverage_estimator}
\title{Filter samples based on ecological coverage estimators}
\description{
Filters samples based on Good's or Chao's coverage estimates, which measure
the completeness of sampling rather than just absolute sequencing depth.
This approach removes samples that are likely undersampled ecologically.
}
\usage{
filter_by_coverage_estimator(table, method = c("good", "chao"),
                              target_coverage = NULL, min_reads = 0)
}
\arguments{
  \item{table}{A feature table (data.frame or matrix) with features as rows and samples as columns.
              First column should be feature IDs, remaining columns are sample counts.}
  \item{method}{Coverage estimation method: "good" for Good's coverage, "chao" for Chao's coverage.}
  \item{target_coverage}{Target coverage threshold (0-1). Samples below this will be filtered.
                        Default is 0.95 for Good's, 0.90 for Chao's.}
  \item{min_reads}{Optional minimum absolute read count cutoff applied in addition to coverage.
                  Useful as a safety floor. Default is 0 (no additional cutoff).}
}
\value{
A list containing:
  \item{table}{Filtered feature table}
  \item{coverage_before}{Coverage estimates before filtering}
  \item{coverage_after}{Coverage estimates after filtering}
  \item{mean_coverage_before}{Mean coverage before filtering}
  \item{mean_coverage_after}{Mean coverage after filtering}
  \item{n_samples_before}{Number of samples before filtering}
  \item{n_samples_after}{Number of samples after filtering}
  \item{n_samples_filtered}{Number of samples removed}
  \item{method}{Method used ("good" or "chao")}
  \item{target_coverage}{Target coverage threshold used}
}
\examples{
# Filter using Good's coverage (keep samples with >= 95\% coverage)
# result <- filter_by_coverage_estimator(table, method = "good", target_coverage = 0.95)

# Filter using Chao's coverage (keep samples with >= 90\% coverage)
# result <- filter_by_coverage_estimator(table, method = "chao", target_coverage = 0.90)
}
