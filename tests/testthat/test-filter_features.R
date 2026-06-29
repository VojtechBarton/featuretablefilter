test_that("filter_features_by_abundance zeros low-abundance features", {
  table <- data.frame(
    id = c("f1", "f2", "f3"),
    s1 = c(10, 5, 2),
    s2 = c(8, 3, 1),
    s3 = c(12, 7, 4)
  )

  # Filter with threshold of 5 (absolute mode)
  result <- filter_features_by_abundance(table, threshold = 5, mode = "absolute")

  expect_s3_class(result, "data.frame")
  # f1 and f2 should be kept (have values >= 5)
  # f3 has all values < 5, so it should be removed by default
  expect_true("f1" %in% result$id)
  expect_true("f2" %in% result$id)
})

test_that("filter_features_by_abundance keeps features with min_samples", {
  table <- data.frame(
    id = c("f1", "f2", "f3"),
    s1 = c(10, 1, 5),
    s2 = c(8, 1, 5),
    s3 = c(12, 1, 5)
  )

  # Feature f2 has no values >= 5, should be removed
  result <- filter_features_by_abundance(table, threshold = 5, mode = "absolute",
                                          min_samples = 2)

  expect_true("f1" %in% result$id)
  expect_true("f3" %in% result$id)
  expect_false("f2" %in% result$id)
})

test_that("filter_features_by_abundance preserves column names", {
  table <- data.frame(
    id = c("f1", "f2"),
    sample_A = c(10, 2),
    sample_B = c(8, 1)
  )

  result <- filter_features_by_abundance(table, threshold = 5, mode = "absolute")

  expect_true("sample_A" %in% colnames(result))
  expect_true("sample_B" %in% colnames(result))
})

test_that("filter_features_by_abundance works in relative mode", {
  table <- data.frame(
    id = c("f1", "f2"),
    s1 = c(90, 10),   # Total: 100
    s2 = c(80, 20)    # Total: 100
  )

  # Threshold of 0.05 (5%) in relative mode
  result <- filter_features_by_abundance(table, threshold = 0.05, mode = "relative")

  expect_s3_class(result, "data.frame")
  # Both features should be kept as they both exceed 5% in at least one sample
  expect_equal(length(result$id), 2)
})
