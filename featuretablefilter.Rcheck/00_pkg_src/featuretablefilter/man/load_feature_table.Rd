\name{load_feature_table}
\alias{load_feature_table}
\title{Load a feature table from TSV or CSV file}
\description{
Automatically detects the file format and structure of microbiome feature tables.
Handles common formats where the first column contains feature IDs.
}
\usage{
load_feature_table(file, header = NULL, sep = NULL,
                    feature_col = 1, feature_name_pattern = NULL)
}
\arguments{
  \item{file}{Path to the input file (TSV or CSV).}
  \item{header}{Logical. Should the first row be treated as headers (sample names)?
              If NULL, will be auto-detected.}
  \item{sep}{File separator. If NULL, will be auto-detected (.tsv uses tab, .csv uses comma).}
  \item{feature_col}{Name or index of the column containing feature IDs. Default is 1 (first column).}
  \item{feature_name_pattern}{Optional regex pattern to identify the feature ID column by name.
                             Common patterns: "OTU", "feature", "ASV", "taxon".}
}
\value{
A data.frame with the feature table. Feature IDs are kept as a regular column
(not row names) to preserve the standard OTU/feature table format.
}
\examples{
# Load TSV file (auto-detects format)
# table <- load_feature_table("feature-table.tsv")

# Load CSV file
# table <- load_feature_table("feature-table.csv")

# Specify feature column by name
# table <- load_feature_table("table.tsv", feature_col = "OTU_ID")
}
