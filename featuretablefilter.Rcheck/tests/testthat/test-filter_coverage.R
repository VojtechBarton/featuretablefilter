test_that("filter_by_coverage removes low-coverage samples", {
  # Create test data
  table <- data.frame(
    id = c("f1", "f2", "f3"),
    s1 = c(100, 50, 30),   # Total: 180
    s2 = c(10, 5, 5),      # Total: 20
    s3 = c(200, 100, 50)   # Total: 350
  )

  # Filter with threshold of 100 reads
  result <- filter_by_coverage(table, min_reads = 100)

  expect_s3_class(result, "data.frame")
  # s2 should be removed (only 20 reads)
  expect_true("s1" %in% colnames(result))
  expect_true("s3" %in% colnames(result))
  expect_false("s2" %in% colnames(result))
  expect_equal(ncol(result), 3)  # id + s1 + s3
})

test_that("filter_by_coverage keeps all samples when threshold is low", {
  table <- data.frame(
    id = c("f1", "f2"),
    s1 = c(100, 50),
    s2 = c(80, 40)
  )

  result <- filter_by_coverage(table, min_reads = 10)

  expect_equal(ncol(result), 3)  # All samples kept
})

test_that("filter_by_coverage preserves column names", {
  table <- data.frame(
    id = c("f1", "f2"),
    sample_A = c(100, 50),
    sample_B = c(80, 40)
  )

  result <- filter_by_coverage(table, min_reads = 10)

  expect_true("sample_A" %in% colnames(result))
  expect_true("sample_B" %in% colnames(result))
})
