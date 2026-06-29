\name{convert_to_relative}
\alias{convert_to_relative}
\title{Convert absolute abundance table to relative abundance table}
\description{
Transforms a feature table from absolute counts (reads) to relative abundances
(proportions) by normalizing each sample to sum to 1 (or 100 if percent = TRUE).
}
\usage{
convert_to_relative(table, percent = FALSE, zero_replace = 0)
}
\arguments{
  \item{table}{A feature table (data.frame or matrix) with features as rows and samples as columns.
              First column should be feature IDs, remaining columns are sample counts.}
  \item{percent}{Logical. If TRUE, return percentages (sum to 100). If FALSE (default),
                return proportions (sum to 1).}
  \item{zero_replace}{Value to replace NA/NaN values that may arise from dividing 0/0.
                     Default is 0.}
}
\value{
A data.frame with the same structure as input, but with relative abundances
instead of absolute counts.
}
\examples{
# Convert to proportions
# rel_table <- convert_to_relative(table)

# Convert to percentages
# rel_table <- convert_to_relative(table, percent = TRUE)
}
