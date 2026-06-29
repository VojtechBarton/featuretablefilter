\name{plot_reads_vs_asvs}
\alias{plot_reads_vs_asvs}
\title{Plot Total Reads vs Observed ASVs with outlier detection}
\description{
Simple visualization showing the relationship between sequencing depth (total reads)
and observed richness (number of non-zero ASVs/features). Uses MAD-based outlier
detection to flag samples that have unusually low richness for their sequencing depth.

This plot helps identify:
\itemize{
  \item Low-depth samples with poor feature discovery (potential undersampling)
  \item High-depth samples with unexpectedly few features (potential technical artifacts)
  \item Samples that deviate from the expected depth-richness curve
}
}
\usage{
plot_reads_vs_asvs(table, mad_multiplier = 3,
                    main = "Total Reads vs Observed ASVs",
                    color = "blue", show_labels = NULL)
}
\arguments{
  \item{table}{A feature table (data.frame or matrix) with features as rows and samples as columns.
              First column should be feature IDs, remaining columns are sample counts.}
  \item{mad_multiplier}{Multiplier for MAD-based outlier detection. Samples with residuals
                       beyond (median + multiplier * MAD) are flagged. Default is 3.}
  \item{main}{Main title for the plot. Default is "Total Reads vs Observed ASVs".}
  \item{color}{Theme color. Default is "blue". Options: "blue", "green", "red", "purple"}
  \item{show_labels}{Logical. Show sample names for outliers. Default is TRUE if <= 10 outliers.}
}
\value{
A list containing:
  \item{plot}{ggplot object}
  \item{outliers}{data.frame of flagged outlier samples}
  \item{metrics}{data.frame with all sample metrics}
}
\examples{
# Simple plot with automatic outlier detection
# result <- plot_reads_vs_asvs(my_table)
# print(result$plot)

# Adjust sensitivity
# result <- plot_reads_vs_asvs(my_table, mad_multiplier = 2.5)

# Access outlier information
# outlier_samples <- result$outliers$sample_name
}
