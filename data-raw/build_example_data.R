# Build example datasets for the package
# This script is not part of the package build; run manually when data change.

extdata <- system.file("extdata", package = "featuretablefilter")

# Feature table (data.frame)
example_feature_table <- read.table(
  file.path(extdata, "example_feature_table.tsv"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# phyloseq and TreeSummarizedExperiment objects (require respective packages)
if (requireNamespace("phyloseq", quietly = TRUE)) {
  example_phyloseq <- readRDS(file.path(extdata, "example_phyloseq_object.rds"))
  save(example_phyloseq, file = "data/example_phyloseq.rda", compress = "xz")
}

if (requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
  example_treesummarizedexperiment <- readRDS(
    file.path(extdata, "example_treesummarizedexperiment_object.rds")
  )
  save(example_treesummarizedexperiment,
       file = "data/example_treesummarizedexperiment.rda", compress = "xz")
}

save(example_feature_table, file = "data/example_feature_table.rda", compress = "xz")
