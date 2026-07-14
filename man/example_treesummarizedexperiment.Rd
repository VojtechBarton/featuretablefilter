\name{example_treesummarizedexperiment}
\alias{example_treesummarizedexperiment}
\docType{data}
\title{Example TreeSummarizedExperiment Object for Testing}

\description{
A synthetic TreeSummarizedExperiment object containing count data, rowData, colData, and a phylogenetic tree.
This dataset is provided for testing the \code{featuretablefilter} package functions
that work with TreeSummarizedExperiment objects.
}

\format{
A TreeSummarizedExperiment object with:
\itemize{
  \item \code{assays}: List containing "counts" matrix (100 features × 30 samples)
  \item \code{rowData}: Feature metadata (Kingdom, Phylum, Class, Mean_Abundance, Prevalence)
  \item \code{colData}: Sample metadata (Environment, pH, Temperature, Treatment)
  \item \code{rowTree}: Phylogenetic tree with 100 tips
  \item \code{rowLink}: Links between features and tree tips
}
}

\source{
Synthetic data generated for package testing and demonstration purposes.
}

\usage{
data(example_treesummarizedexperiment)
}

\examples{
# Load the example TreeSummarizedExperiment object
data(example_treesummarizedexperiment)

# View basic structure
tse <- example_treesummarizedexperiment
print(tse)

# Run pipeline on TSE object (returns TSE)
library(featuretablefilter)
result <- run_filtering_pipeline(
  input = tse,
  prefix = "filtered",
  cov_filter_method = "iqr",
  cov_threshold = 1.5,
  abun_filter_method = "relative_cutoff",
  abun_threshold = 0.01
)

# Result is still a TreeSummarizedExperiment
filtered_tse <- result$filtered_table
class(filtered_tse)  # "TreeSummarizedExperiment"
}

\seealso{
\code{\link{run_filtering_pipeline}},
\code{\link{from_TSE}},
\code{\link{to_TSE}}
}

\keyword{datasets}
