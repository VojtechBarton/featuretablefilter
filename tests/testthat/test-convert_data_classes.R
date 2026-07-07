test_that("from_phyloseq converts correctly", {
  skip_if_not_installed("phyloseq")

  library(phyloseq)

  # Create test phyloseq object
  otu_mat <- matrix(rpois(50, 5), nrow = 10, ncol = 5)
  rownames(otu_mat) <- paste0("ASV_", 1:10)
  colnames(otu_mat) <- paste0("Sample_", 1:5)

  ps <- phyloseq(
    otu_table(otu_mat, taxa_are_rows = TRUE)
  )

  # Convert to data.frame
  df <- from_phyloseq(ps)

  expect_s3_class(df, "data.frame")
  expect_true("feature_id" %in% colnames(df))
  expect_equal(nrow(df), 10)
  expect_equal(ncol(df), 6)  # feature_id + 5 samples
})

test_that("from_phyloseq with taxonomy", {
  skip_if_not_installed("phyloseq")

  library(phyloseq)

  # Create test phyloseq with taxonomy
  otu_mat <- matrix(rpois(30, 5), nrow = 5, ncol = 6)
  rownames(otu_mat) <- paste0("ASV_", 1:5)
  colnames(otu_mat) <- paste0("Sample_", 1:6)

  tax_mat <- matrix(c(
    "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Genus1",
    "Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Bacteroidaceae", "Genus2"
  ), nrow = 5, ncol = 6, byrow = TRUE)
  rownames(tax_mat) <- paste0("ASV_", 1:5)
  colnames(tax_mat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

  ps <- phyloseq(
    otu_table(otu_mat, taxa_are_rows = TRUE),
    tax_table(tax_mat)
  )

  # Convert with taxonomy
  df <- from_phyloseq(ps, include_taxa = TRUE)

  expect_s3_class(df, "data.frame")
  expect_true("feature_id" %in% colnames(df))
  expect_true("tax_Kingdom" %in% colnames(df))
  expect_true("tax_Genus" %in% colnames(df))
})

test_that("to_phyloseq creates valid object", {
  skip_if_not_installed("phyloseq")

  library(phyloseq)

  # Create test data.frame
  df <- data.frame(
    feature_id = paste0("ASV_", 1:10),
    Sample_1 = rpois(10, 5),
    Sample_2 = rpois(10, 5),
    Sample_3 = rpois(10, 5)
  )

  # Convert to phyloseq
  ps <- to_phyloseq(df)

  expect_s4_class(ps, "phyloseq")
  expect_true(!is.null(otu_table(ps)))
  expect_true(taxa_are_rows(otu_table(ps)))
  expect_equal(nrow(otu_table(ps)), 10)
  expect_equal(ncol(otu_table(ps)), 3)
})

test_that("to_phyloseq with taxonomy", {
  skip_if_not_installed("phyloseq")

  library(phyloseq)

  # Create test data
  df <- data.frame(
    feature_id = paste0("ASV_", 1:5),
    S1 = rpois(5, 5),
    S2 = rpois(5, 5)
  )

  tax_df <- data.frame(
    feature_id = paste0("ASV_", 1:5),
    Kingdom = rep("Bacteria", 5),
    Phylum = c("Firmicutes", "Bacteroidetes", "Proteobacteria", "Actinobacteria", "Firmicutes")
  )

  ps <- to_phyloseq(df, tax_table = tax_df)

  expect_s4_class(ps, "phyloseq")
  expect_true(!is.null(tax_table(ps)))
  expect_equal(nrow(tax_table(ps)), 5)
})

test_that("from_TSE converts correctly", {
  skip_if_not_installed("TreeSummarizedExperiment")

  library(TreeSummarizedExperiment)
  library(SummarizedExperiment)

  # Create test TSE
  assay_mat <- matrix(rpois(50, 5), nrow = 10, ncol = 5)
  rownames(assay_mat) <- paste0("ASV_", 1:10)
  colnames(assay_mat) <- paste0("Sample_", 1:5)

  tse <- TreeSummarizedExperiment(
    assays = list(counts = assay_mat),
    rowData = DataFrame(Genus = paste0("Genus_", 1:10)),
    colData = DataFrame(Condition = rep(c("A", "B"), c(3, 2)))
  )

  # Convert to data.frame
  df <- from_TSE(tse)

  expect_s3_class(df, "data.frame")
  expect_true("feature_id" %in% colnames(df))
  expect_equal(nrow(df), 10)
})

test_that("from_TSE with rowData and colData", {
  skip_if_not_installed("TreeSummarizedExperiment")

  library(TreeSummarizedExperiment)

  assay_mat <- matrix(rpois(30, 5), nrow = 6, ncol = 5)
  rownames(assay_mat) <- paste0("ASV_", 1:6)
  colnames(assay_mat) <- paste0("Sample_", 1:5)

  tse <- TreeSummarizedExperiment(
    assays = list(counts = assay_mat),
    rowData = DataFrame(Kingdom = rep("Bacteria", 6), Phylum = paste0("Phylum_", 1:6)),
    colData = DataFrame(Condition = rep(c("Control", "Treatment"), c(2, 3)))
  )

  # Convert with rowData
  df <- from_TSE(tse, add_row_data = TRUE)

  expect_true("rowData_Kingdom" %in% colnames(df))
  expect_true("rowData_Phylum" %in% colnames(df))

  # Convert with colData
  result <- from_TSE(tse, add_col_data = TRUE)

  expect_type(result, "list")
  expect_true(all(c("table", "sample_data") %in% names(result)))
  expect_s3_class(result$sample_data, "data.frame")
})

test_that("to_TSE creates valid object", {
  skip_if_not_installed("TreeSummarizedExperiment")

  library(TreeSummarizedExperiment)

  # Create test data.frame
  df <- data.frame(
    feature_id = paste0("ASV_", 1:10),
    Sample_1 = rpois(10, 5),
    Sample_2 = rpois(10, 5),
    Sample_3 = rpois(10, 5)
  )

  # Convert to TSE
  tse <- to_TSE(df)

  expect_s4_class(tse, "TreeSummarizedExperiment")
  expect_true("counts" %in% names(assays(tse)))
  expect_equal(nrow(tse), 10)
  expect_equal(ncol(tse), 3)
})

test_that("to_TSE with rowData and colData", {
  skip_if_not_installed("TreeSummarizedExperiment")

  library(TreeSummarizedExperiment)

  df <- data.frame(
    feature_id = paste0("ASV_", 1:6),
    S1 = rpois(6, 5),
    S2 = rpois(6, 5),
    S3 = rpois(6, 5)
  )

  rowData <- data.frame(
    Kingdom = rep("Bacteria", 6),
    Phylum = paste0("Phylum_", 1:6)
  )
  rownames(rowData) <- paste0("ASV_", 1:6)

  colData <- data.frame(
    Condition = rep(c("Control", "Treatment"), 3)
  )
  rownames(colData) <- c("S1", "S2", "S3")

  tse <- to_TSE(df, rowData = rowData, colData = colData)

  expect_true(!is.null(rowData(tse)))
  expect_true(!is.null(colData(tse)))
  expect_equal(nrow(rowData(tse)), 6)
  expect_equal(nrow(colData(tse)), 3)
})

test_that("convert_feature_table works for all conversions", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("TreeSummarizedExperiment")

  library(phyloseq)
  library(TreeSummarizedExperiment)

  # Start with data.frame
  df <- data.frame(
    feature_id = paste0("ASV_", 1:10),
    S1 = rpois(10, 5),
    S2 = rpois(10, 5)
  )

  # Convert to phyloseq
  ps <- convert_feature_table(df, to = "phyloseq")
  expect_s4_class(ps, "phyloseq")

  # Convert back to data.frame
  df_back <- convert_feature_table(ps, to = "data.frame")
  expect_s3_class(df_back, "data.frame")
  expect_equal(nrow(df_back), 10)

  # Convert to TSE
  tse <- convert_feature_table(df, to = "TSE")
  expect_s4_class(tse, "TreeSummarizedExperiment")

  # Convert TSE to data.frame
  df_from_tse <- convert_feature_table(tse, to = "data.frame")
  expect_s3_class(df_from_tse, "data.frame")
})

test_that("convert_feature_table handles unsupported classes", {
  expect_error(
    convert_feature_table(list(a = 1), to = "data.frame"),
    "Cannot convert"
  )
})

test_that("from_phyloseq handles OTU table orientation", {
  skip_if_not_installed("phyloseq")

  library(phyloseq)

  # Create phyloseq with taxa_are_rows = FALSE
  otu_mat <- matrix(rpois(30, 5), nrow = 6, ncol = 5)
  rownames(otu_mat) <- paste0("Sample_", 1:6)
  colnames(otu_mat) <- paste0("ASV_", 1:5)

  ps <- phyloseq(
    otu_table(otu_mat, taxa_are_rows = FALSE)
  )

  df <- from_phyloseq(ps)

  expect_equal(nrow(df), 5)  # Features should be rows
  expect_true(all(grepl("ASV", df$feature_id)))
})

test_that("to_phyloseq validates input", {
  skip_if_not_installed("phyloseq")

  expect_error(
    to_phyloseq(matrix(1:10)),
    "Could not detect feature ID"
  )
})

test_that("from_TSE validates input", {
  skip_if_not_installed("TreeSummarizedExperiment")

  expect_error(
    from_TSE(list(a = 1)),
    "Input must be a"
  )
})

test_that("to_TSE with custom assay name", {
  skip_if_not_installed("TreeSummarizedExperiment")

  df <- data.frame(
    feature_id = paste0("ASV_", 1:5),
    S1 = rpois(5, 5),
    S2 = rpois(5, 5)
  )

  tse <- to_TSE(df, assay_name = "abundance")

  expect_true("abundance" %in% names(assays(tse)))
  expect_false("counts" %in% names(assays(tse)))
})
