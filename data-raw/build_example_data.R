# Build example datasets for the package
# This script is not part of the package build; run manually when data change.

extdata <- system.file("extdata", package = "featuretablefilter")

# Feature table (data.frame) — saved to data/ for lazy loading via data()
example_feature_table <- read.table(
  file.path(extdata, "example_feature_table.tsv"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

save(example_feature_table, file = "data/example_feature_table.rda", compress = "xz")

# Note: example_phyloseq_object.rds and example_treesummarizedexperiment_object.rds
# are stored in inst/extdata/ as .rds files because they require the respective
# packages to load. They are accessed via system.file() at runtime, not via data().
