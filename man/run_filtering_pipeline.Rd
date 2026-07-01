\name{run_filtering_pipeline}
\alias{run_filtering_pipeline}
\title{Run a complete feature table filtering pipeline with QC report}
\description{
A comprehensive pipeline that loads, filters, and analyzes a feature table.
Allows selection of different coverage and abundance filtering methods.
Produces QC metrics and visualizations comparing before/after filtering.

Generates the following outputs: filtered feature table (TSV format), sample coverage distribution plot, feature abundance distribution plot, sample sparsity histogram (per-sample view), retention rates barplot, heatmaps of top variable features (original and filtered), top features stacked barplot, presence frequency comparison plot (side-by-side histograms), text summary report.
}
\usage{
run_filtering_pipeline(input_file,
                        output_dir = ".",
                        prefix = "filtered",
                        # Coverage filtering
                        cov_filter_method = c("none", "absolute", "mad", "iqr"),
                        cov_threshold = NULL,
                        cov_floor = 0,
                        # Abundance filtering
                        abun_filter_method = c("none", "absolute", "relative", "relative_cutoff"),
                        abun_threshold = NULL,
                        abun_min_samples = 1,
                        # Additional options
                        min_coverage_for_relative = 1000,
                        remove_features = TRUE,
                        generate_plots = TRUE,
                        generate_report = TRUE,
                        verbose = TRUE)
}
\arguments{
  \item{input_file}{Path to input feature table file (TSV or CSV).}
  \item{output_dir}{Directory to save outputs (filtered table, plots, report). Default is current directory.}
  \item{prefix}{Prefix for output filenames. Default is "filtered".}
  \item{cov_filter_method}{Method for sample coverage filtering:
                          \describe{
                            \item{"none"}{No coverage filtering}
                            \item{"absolute"}{Use fixed min_reads threshold}
                            \item{"mad"}{Use MAD-based threshold estimation (median +/- multiplier * MAD)}
                            \item{"iqr"}{Use IQR/Tukey-based threshold estimation (fences)}
                          }}
  \item{cov_threshold}{Fixed threshold if method="absolute", or multiplier if method="mad"/"iqr".
                      For MAD/IQR, this multiplies the spread statistic.}
  \item{cov_floor}{Minimum possible coverage cutoff when using MAD/IQR methods.}
  \item{abun_filter_method}{Method for feature abundance filtering:
                           \describe{
                             \item{"none"}{No abundance filtering}
                             \item{"absolute"}{Filter features below absolute count threshold}
                             \item{"relative"}{Filter features below relative abundance threshold}
                             \item{"relative_cutoff"}{Use relative threshold based on min-coverage sample}
                           }}
  \item{abun_threshold}{Threshold value for abundance filtering.
                       For "absolute": minimum read count per feature.
                       For "relative": proportion (e.g., 0.01 for 1\%).
                       For "relative_cutoff": relative proportion used to calculate absolute threshold.}
  \item{abun_min_samples}{Minimum number of samples where feature must exceed threshold (default: 1).}
  \item{min_coverage_for_relative}{If using "relative_cutoff" method, minimum coverage for samples.}
  \item{remove_features}{Logical. Remove features below threshold? Default TRUE.}
  \item{generate_plots}{Logical. Generate QC visualization plots? Default TRUE.}
  \item{generate_report}{Logical. Generate text summary report? Default TRUE.}
  \item{verbose}{Logical. Print progress messages? Default TRUE.}
}
\value{
A list containing:
  \item{original_table}{Original loaded table}
  \item{filtered_table}{Filtered table}
  \item{qc_metrics}{List of QC metrics from compute_filtering_qc()}
  \item{presence_stats}{List of presence analysis statistics (including comparison data)}
  \item{filtering_summary}{Summary of filtering steps applied}
}
\examples{
# Simple pipeline with absolute thresholds
# result <- run_filtering_pipeline(
#   input_file = "feature-table.tsv",
#   output_dir = "results",
#   prefix = "analysis1",
#   cov_filter_method = "absolute",
#   cov_threshold = 1000,
#   abun_filter_method = "absolute",
#   abun_threshold = 5
# )

# Pipeline with MAD-based coverage and relative abundance filtering
# result <- run_filtering_pipeline(
#   input_file = "feature-table.tsv",
#   output_dir = "results",
#   cov_filter_method = "mad",
#   cov_threshold = 3,
#   abun_filter_method = "relative_cutoff",
#   abun_threshold = 0.01,
#   min_coverage_for_relative = 1000
# )
}
