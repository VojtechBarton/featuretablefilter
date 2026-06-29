\name{plot_coverage_histogram}
\alias{plot_coverage_histogram}
\title{Plot histogram of sample coverage}
\description{
Creates a histogram of total reads per sample with dynamic bin sizing.
Optionally highlights samples above/below a threshold.
}
\usage{
plot_coverage_histogram(table, threshold = NULL,
                         main = "Sample Coverage Distribution",
                         xlab = "Total Reads per Sample",
                         col_above = "steelblue",
                         col_below = "tomato",
                         plot_dir = NULL, prefix = "coverage_histogram",
                         width = 8, height = 6, dpi = 300)
}
\arguments{
  \item{table}{A feature table (data.frame or matrix) with features as rows and samples as columns.
              First column should be feature IDs, remaining columns are sample counts.}
  \item{threshold}{Optional coverage threshold. Bars below this value are colored red (to remove),
                  bars at or above are colored green (to keep).}
  \item{main}{Title for the plot. Default is "Sample Coverage Distribution".}
  \item{xlab}{Label for x-axis. Default is "Total Reads per Sample".}
  \item{col_above}{Color for bars at or above threshold (default: "steelblue").}
  \item{col_below}{Color for bars below threshold (default: "tomato").}
  \item{plot_dir}{Optional directory path to save the plot. If NULL, plot is displayed only.}
  \item{prefix}{Prefix for saved plot filename. Default is "coverage_histogram".}
  \item{width}{Plot width in inches for ggsave. Default 8.}
  \item{height}{Plot height in inches for ggsave. Default 6.}
  \item{dpi}{DPI for saved raster output. Default 300.}
}
\value{
Invisibly returns a list with:
  \item{coverage}{Vector of total reads per sample}
  \item{breaks}{The breaks used for the histogram}
  \item{counts}{Counts in each bin}
  \item{plot_path}{Path to saved plot (if plot_dir provided)}
  \item{plot}{The ggplot object}
}
\examples{
# Basic histogram
# plot_coverage_histogram(table)

# With threshold highlighting
# plot_coverage_histogram(table, threshold = 1000)

# Save plot to file
# plot_coverage_histogram(table, threshold = 1000, plot_dir = "plots", prefix = "my_plot")
}
