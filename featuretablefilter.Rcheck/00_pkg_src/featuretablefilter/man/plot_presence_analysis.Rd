\name{plot_presence_analysis}
\alias{plot_presence_analysis}
\title{Plot presence frequency analysis for features and samples}
\description{
Creates histograms showing:
\itemize{
  \item In how many samples each feature is present (feature prevalence)
  \item How many features are present in each sample (sample richness)
}

When \code{table_filtered} is provided, generates side-by-side comparison plots with shared axes.
}
\usage{
plot_presence_analysis(table, threshold = 1,
                        save_dir = NULL, prefix = "presence_analysis",
                        main_main = "Presence Frequency Analysis",
                        width = 10, height = 6, dpi = 300,
                        prev_bins = 20, rich_bins = 20,
                        table_filtered = NULL, prefix_filtered = "filtered")
}
\arguments{
  \item{table}{A feature table (data.frame or matrix) with features as rows and samples as columns.
              First column should be feature IDs, remaining columns are sample counts.}
  \item{threshold}{Minimum abundance for a feature to be considered "present" in a sample.
                  Default is 1 (any non-zero count).}
  \item{save_dir}{Optional directory to save plots. If NULL, plots are displayed only.}
  \item{prefix}{Prefix for saved plot filenames. Default is "presence_analysis".}
  \item{main_main}{Main title for the combined plot. Default is "Presence Frequency Analysis".}
  \item{width}{Plot width in inches for ggsave. Default 10.}
  \item{height}{Plot height in inches for ggsave. Default 6.}
  \item{dpi}{DPI for saved raster output. Default 300.}
  \item{prev_bins}{Number of bins for feature prevalence histogram. Default is 20.}
  \item{rich_bins}{Number of bins for sample richness histogram. Default is 20.}
  \item{table_filtered}{Optional filtered feature table for comparison. If provided,
                       shows original vs filtered side by side with shared axes.}
  \item{prefix_filtered}{Prefix for filtered plot files. Default is "filtered".}
}
\value{
A list containing:
  \item{feature_prevalence}{Named vector: number of samples each feature appears in (original)}
  \item{sample_richness}{Named vector: number of features in each sample (original)}
  \item{feature_prevalence_filtered}{Named vector: filtered feature prevalence (if table_filtered provided)}
  \item{sample_richness_filtered}{Named vector: filtered sample richness (if table_filtered provided)}
  \item{mean_feature_prevalence}{Mean number of samples per feature (original)}
  \item{mean_sample_richness}{Mean number of features per sample (original)}
  \item{mean_feature_prevalence_filtered}{Mean for filtered (if table_filtered provided)}
  \item{mean_sample_richness_filtered}{Mean for filtered (if table_filtered provided)}
  \item{min_feature_prevalence}{Minimum feature prevalence (original)}
  \item{max_feature_prevalence}{Maximum feature prevalence (original)}
  \item{min_sample_richness}{Minimum sample richness (original)}
  \item{max_sample_richness}{Maximum sample richness (original)}
  \item{min_feature_prevalence_filtered}{Minimum for filtered (if table_filtered provided)}
  \item{max_feature_prevalence_filtered}{Maximum for filtered (if table_filtered provided)}
  \item{min_sample_richness_filtered}{Minimum for filtered (if table_filtered provided)}
  \item{max_sample_richness_filtered}{Maximum for filtered (if table_filtered provided)}
  \item{n_features}{Total number of features (original)}
  \item{n_samples}{Total number of samples (original)}
  \item{n_features_filtered}{Total for filtered (if table_filtered provided)}
  \item{n_samples_filtered}{Total for filtered (if table_filtered provided)}
  \item{plots}{List containing the ggplot objects}
}
\examples{
# Display plots for original table only
# plot_presence_analysis(table)

# Compare original vs filtered (side-by-side with shared axes)
# plot_presence_analysis(table, table_filtered = filtered_table)

# Save plots
# plot_presence_analysis(table, save_dir = "plots", prefix = "my_analysis")

# Get statistics
# stats <- plot_presence_analysis(table)
# print(stats$mean_sample_richness)
}
