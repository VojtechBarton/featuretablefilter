\name{plot_qc_comparison}
\alias{plot_qc_comparison}
\title{Create visualization comparison of original vs filtered feature tables}
\description{
Generates a comprehensive set of plots to compare original and filtered tables.
Includes coverage distributions, sparsity histograms, retention barplots, top features barplots, and heatmaps.

Plots generated:
\itemize{
  \item Sample coverage distribution (log10 scale, faceted)
  \item Feature abundance distribution (log10 scale, faceted)
  \item Sample sparsity histogram (per-sample sparsity distribution, faceted with mean lines)
  \item Retention rates barplot (read, feature, sample retention)
  \item Top N features stacked barplot (relative abundance comparison)
  \item Heatmap of top variable features (original and filtered separately)
}
}
\usage{
plot_qc_comparison(original_table, filtered_table,
                    plot_dir = NULL, prefix = "qc_comparison",
                    top_n = 10, heatmap_top = 50,
                    width = 10, height = 8, dpi = 300)
}
\arguments{
  \item{original_table}{Original feature table (data.frame) before filtering.}
  \item{filtered_table}{Filtered feature table (data.frame) after filtering.}
  \item{plot_dir}{Optional directory path to save plots. If NULL, plots are displayed only.}
  \item{prefix}{Prefix for saved plot filenames. Default is "qc_comparison".}
  \item{top_n}{Number of top features to show in barplots. Default is 10.}
  \item{heatmap_top}{Number of top variable features to show in heatmap. Default is 50.}
  \item{width}{Plot width in inches for ggsave. Default 10.}
  \item{height}{Plot height in inches for ggsave. Default 8.}
  \item{dpi}{DPI for saved raster output. Default 300.}
}
\value{
A list containing:
\item{plots}{Named list of ggplot objects}
\item{paths}{Named list of file paths (if plot_dir provided)}
}
\seealso{
\code{\link{compute_filtering_qc}} for numerical QC metrics.
\code{\link{plot_presence_analysis}} for presence frequency analysis.
}
\examples{
# Display plots
# plot_qc_comparison(original_table, filtered_table)

# Save plots to directory
# plot_qc_comparison(original_table, filtered_table, plot_dir = "qc_plots", prefix = "my_analysis")
}
