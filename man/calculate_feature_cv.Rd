\name{calculate_feature_cv}
\alias{calculate_feature_cv}
\title{Calculate coefficient of variation (CV) for features}
\description{
Computes the coefficient of variation (standard deviation / mean) for each feature
across samples. This is a useful metric for identifying highly variable features.
}
\usage{
calculate_feature_cv(table, plot = FALSE, plot_type = c("histogram", "boxplot"),
                      log_transform = FALSE, main = "Feature Coefficient of Variation")
}
\arguments{
  \item{table}{A feature table (data.frame or matrix) with features as rows and samples as columns.
              First column should be feature IDs, remaining columns are sample counts/abundances.}
  \item{plot}{Logical. If TRUE, creates a histogram of CV distribution. Default is FALSE.}
  \item{plot_type}{Type of plot: "histogram" (default) or "boxplot".}
  \item{log_transform}{Logical. If TRUE, applies log10 transformation to CV values for plotting.
                      Useful when CV values span several orders of magnitude. Default is FALSE.}
  \item{main}{Plot title (only used if plot = TRUE).}
}
\value{
A data.frame with feature IDs and their corresponding CV values.
        Columns: feature_id, mean, sd, cv
}
\examples{
# Calculate CV for all features
# cv_results <- calculate_feature_cv(table)

# Calculate and plot distribution
# cv_results <- calculate_feature_cv(table, plot = TRUE)

# Find most variable features
# cv_results <- calculate_feature_cv(table)
# top_variable <- cv_results[order(-cv_results$cv), ][1:10, ]
}
