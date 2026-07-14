\name{example_phyloseq}
\alias{example_phyloseq}
\docType{data}
\title{Example phyloseq Object for Testing}

\description{
A synthetic phyloseq object containing OTU counts, taxonomy, sample data, and a phylogenetic tree.
This dataset is provided for testing the \code{featuretablefilter} package functions
that work with phyloseq objects.
}

\format{
A phyloseq object with:
\itemize{
  \item \code{otu_table}: 100 features (ASVs) × 30 samples
  \item \code{tax_table}: Taxonomy assignments (Kingdom through Species)
  \item \code{sample_data}: Sample metadata (Environment, pH, Temperature, Treatment)
  \item \code{phy_tree}: Random phylogenetic tree with 100 tips
}
}

\source{
Synthetic data generated for package testing and demonstration purposes.
}

\usage{
data(example_phyloseq)
}

\examples{
# Load the example phyloseq object
data(example_phyloseq)

# View basic structure
ps <- example_phyloseq
print(ps)

# Run pipeline on phyloseq object (returns phyloseq)
library(featuretablefilter)
result <- run_filtering_pipeline(
  input = ps,
  prefix = "filtered",
  cov_filter_method = "mad",
  cov_threshold = 2,
  abun_filter_method = "relative_cutoff",
  abun_threshold = 0.01
)

# Result is still a phyloseq object
filtered_ps <- result$filtered_table
class(filtered_ps)  # "phyloseq"
}

\seealso{
\code{\link{run_filtering_pipeline}},
\code{\link{from_phyloseq}},
\code{\link{to_phyloseq}}
}

\keyword{datasets}
