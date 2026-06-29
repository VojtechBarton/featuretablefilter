\name{filter_by_coverage}
\alias{filter_by_coverage}
\title{Filter samples by minimum coverage}
\description{
Removes samples with total read count below the threshold.
}
\usage{
filter_by_coverage(table, min_reads)
}
\arguments{
  \item{table}{A feature table (data.frame or matrix) with features as rows and samples as columns.
              First column should be feature IDs, remaining columns are sample counts.}
  \item{min_reads}{Minimum number of reads required to keep a sample.}
}
\value{
A filtered feature table (data.frame) with low-coverage samples removed.
}
\examples{
# filter_by_coverage(my_table, min_reads = 1000)
}
