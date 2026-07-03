test_that("load_feature_table loads TSV files correctly", {
  # Create a temporary test file
  tmp_file <- tempfile(fileext = ".tsv")
  writeLines(c("#OTU ID\tsample1\tsample2\tsample3",
               "feat1\t10\t20\t30",
               "feat2\t5\t15\t25"),
             tmp_file)

  # Test loading
  result <- load_feature_table(tmp_file)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 4)
  expect_equal(colnames(result)[1], "#OTU ID")
  expect_equal(result$sample1[1], 10)

  unlink(tmp_file)
})

test_that("load_feature_table handles header detection", {
  # File with header starting with #
  tmp_file <- tempfile(fileext = ".tsv")
  writeLines(c("#Feature ID\ts1\ts2",
               "a\t1\t2",
               "b\t3\t4"),
             tmp_file)

  result <- load_feature_table(tmp_file)

  expect_true("#Feature ID" %in% colnames(result))
  expect_equal(nrow(result), 2)

  unlink(tmp_file)
})
