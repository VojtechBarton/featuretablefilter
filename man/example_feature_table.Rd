\name{example_feature_table}
\alias{example_feature_table}
\docType{data}
\title{Example Feature Table for Testing and Documentation}

\description{
A synthetic microbiome feature table containing 50 features (ASVs) across 20 samples.
This dataset is provided for testing the \code{featuretablefilter} package functions
and for use in documentation examples.
}

\format{
A data frame with 50 rows (features) and 21 columns:
\itemize{
  \item \code{Feature_ID}: Feature identifier (ASV_001 to ASV_050)
  \item \code{Sample_001} to \code{Sample_020}: Read counts for each sample
}

The data simulates typical microbiome sequencing data with:
\itemize{
  \item Variable coverage across samples (ranging from ~5,000 to ~25,000 reads per sample)
  \item Variable abundance across features (from rare to abundant taxa)
  \item Realistic sparsity pattern (many zeros for rare features)
  \item Nested structure where dominant features show higher variance
}
}

\source{
Synthetic data generated for package testing and demonstration purposes.
}

\usage{
data(example_feature_table)
}

\examples{
# Load the example data
data(example_feature_table)

# View basic structure
str(example_feature_table)

# Calculate total reads per sample
sample_totals <- colSums(example_table[, -1])
print(sample_totals)

# Run a simple filtering pipeline
library(featuretablefilter)
result <- run_filtering_pipeline(
  input = example_feature_table,
  cov_filter_method = "mad",
  cov_threshold = 2,
  abun_filter_method = "absolute",
  abun_threshold = 10
)
}

\seealso{
\code{\link{run_filtering_pipeline}},
\code{\link{load_feature_table}}
}

\keyword{datasets}
