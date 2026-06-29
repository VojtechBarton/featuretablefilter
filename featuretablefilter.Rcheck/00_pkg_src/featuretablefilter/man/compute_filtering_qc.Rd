\name{compute_filtering_qc}
\alias{compute_filtering_qc}
\title{Compute quality control metrics comparing original and filtered tables}
\description{
Calculates multiple QC metrics to assess the impact of filtering on a feature table.
Includes sparsity change, retention rates, rank-abundance stability, and compositional similarity.
}
\usage{
compute_filtering_qc(original_table, filtered_table, top_n = 10)
}
\arguments{
  \item{original_table}{Original feature table (data.frame) before filtering.}
  \item{filtered_table}{Filtered feature table (data.frame) after filtering.}
  \item{top_n}{Number of top taxa to consider for Rank-Abundance Stability. Default is 10.}
}
\value{
A list containing:
  \item{sparsity_original}{Sparsity of original table}
  \item{sparsity_filtered}{Sparsity of filtered table}
  \item{sparsity_drop_percent}{Percentage point drop in sparsity}
  \item{read_retention_percent}{Percentage of reads retained}
  \item{feature_retention_percent}{Percentage of features retained}
  \item{sample_retention_percent}{Percentage of samples retained}
  \item{rank_abundance_correlation}{Spearman correlation of top N taxa abundances}
  \item{rank_abundance_pvalue}{P-value for rank-abundance correlation}
  \item{procrustes_m2}{Procrustes M^2 statistic (Bray-Curtis based)}
  \item{procrustes_correlation}{Procrustes correlation coefficient}
}
\examples{
# Compare original and filtered tables
# qc <- compute_filtering_qc(original_table, filtered_table)
# print(qc)
}
