\name{filter_features_joint}
\alias{filter_features_joint}
\title{Filter features using joint abundance and prevalence criteria}
\description{
Combines abundance thresholds with prevalence requirements using AND/OR logic.
This allows flexible filtering where features can be retained based on being
either abundant, widespread, or both.

Common use cases:
\itemize{
  \item AND logic: Keep only features that are BOTH abundant enough AND widely distributed
  \item OR logic: Keep features that are EITHER abundant enough OR widely distributed
    (useful for retaining rare but high-abundance features, or common but low-abundance features)
}
}
\usage{
filter_features_joint(table,
                       abundance_threshold = 0,
                       prevalence_threshold = 0.3,
                       mode = c("relative", "absolute"),
                       logic = c("OR", "AND"),
                       remove_zeros = TRUE)
}
\arguments{
  \item{table}{A feature table (data.frame or matrix) with features as rows and samples as columns.
              First column should be feature IDs, remaining columns are sample counts/abundances.}
  \item{abundance_threshold}{Minimum abundance threshold. Features below this in a sample are considered "absent"
                            for the abundance criterion. Default is 0 (any non-zero value counts).}
  \item{prevalence_threshold}{Proportion of samples (0-1) where feature must meet abundance threshold.
                             For example, 0.3 means the feature must exceed abundance_threshold in at least 30\% of samples.}
  \item{mode}{Abundance measurement mode: "absolute" for raw counts, "relative" for proportions.
             Default is "relative".}
  \item{logic}{Logical operator: "AND" requires both criteria to be met, "OR" requires either criterion.
              Default is "OR".}
  \item{remove_zeros}{Logical. If TRUE (default), remove features that don't meet the filtering criteria.}
}
\value{
A list containing:
  \item{table}{Filtered feature table (data.frame)}
  \item{n_features_before}{Number of features before filtering}
  \item{n_features_after}{Number of features after filtering}
  \item{n_features_removed}{Number of features removed}
  \item{n_by_abundance_only}{Features meeting only abundance criterion (not prevalence)}
  \item{n_by_prevalence_only}{Features meeting only prevalence criterion (not abundance)}
  \item{n_by_both}{Features meeting both criteria}
  \item{n_by_neither}{Features meeting neither criterion (these were removed)}
  \item{abundance_threshold}{The abundance threshold used}
  \item{prevalence_threshold}{The prevalence threshold used}
  \item{logic}{The logic operator used ("AND" or "OR")}
  \item{mode}{The abundance mode used ("absolute" or "relative")}
  \item{feature_details}{Data.frame with per-feature filtering details}
}
\examples{
# Keep features that are EITHER >0.1\% relative abundance OR present in >=30\% of samples
# result <- filter_features_joint(table,
#                                 abundance_threshold = 0.001,
#                                 prevalence_threshold = 0.3,
#                                 mode = "relative",
#                                 logic = "OR")

# Keep features that are BOTH >5 reads AND present in >=50\% of samples
# result <- filter_features_joint(table,
#                                 abundance_threshold = 5,
#                                 prevalence_threshold = 0.5,
#                                 mode = "absolute",
#                                 logic = "AND")
}
