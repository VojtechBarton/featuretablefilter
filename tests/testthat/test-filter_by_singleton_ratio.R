test_that("filter_by_singleton_ratio removes high-ratio samples", {
  # Create a feature table with some suspicious samples
  # Sample A: normal - low singleton ratio
  # Sample B: suspicious - high singleton ratio (many features with count=1)
  # Sample C: normal
  test_table <- data.frame(
    feature_id = c("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8"),
    Sample_A = c(100, 50, 30, 20, 15, 10, 5, 5),   # Total: 235, singletons: 0, ratio: 0
    Sample_B = c(1, 1, 1, 1, 1, 1, 1, 93),         # Total: 100, singletons: 7, ratio: 0.07
    Sample_C = c(80, 40, 20, 10, 5, 3, 2, 1),      # Total: 161, singletons: 1, ratio: 0.006
    stringsAsFactors = FALSE
  )

  # With max_singleton_ratio = 0.05, Sample_B should be removed
  result <- filter_by_singleton_ratio(test_table, max_singleton_ratio = 0.05)

  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 3)  # feature_id + Sample_A + Sample_C
  expect_equal(colnames(result), c("feature_id", "Sample_A", "Sample_C"))
  expect_equal(attr(result, "n_filtered_out"), 1)
  expect_equal(attr(result, "n_retained"), 2)
})

test_that("filter_by_singleton_ratio preserves all samples when below threshold", {
  test_table <- data.frame(
    feature_id = c("F1", "F2", "F3", "F4"),
    Sample_A = c(100, 50, 30, 20),
    Sample_B = c(80, 40, 15, 5),
    stringsAsFactors = FALSE
  )

  # All samples have very low singleton ratios
  result <- filter_by_singleton_ratio(test_table, max_singleton_ratio = 0.1)

  expect_equal(ncol(result), 3)  # feature_id + both samples
  expect_equal(attr(result, "n_filtered_out"), 0)
  expect_equal(attr(result, "n_retained"), 2)
})

test_that("filter_by_singleton_ratio works with count_type='singleton'", {
  test_table <- data.frame(
    feature_id = c("F1", "F2", "F3", "F4", "F5"),
    Sample_A = c(1, 1, 1, 90, 50),   # 3 singletons, 0 doubletons
    Sample_B = c(2, 2, 2, 80, 50),   # 0 singletons, 3 doubletons
    stringsAsFactors = FALSE
  )

  # With count_type="singleton" and threshold 0.05:
  # Sample_A: 3/143 = 0.021 -> kept
  # Sample_B: 0/134 = 0 -> kept
  result_single <- filter_by_singleton_ratio(test_table, max_singleton_ratio = 0.05, count_type = "singleton")
  expect_equal(attr(result_single, "n_retained"), 2)

  # With count_type="doubleton":
  # Sample_A: 0/143 = 0 -> kept
  # Sample_B: 3/134 = 0.022 -> kept
  result_double <- filter_by_singleton_ratio(test_table, max_singleton_ratio = 0.05, count_type = "doubleton")
  expect_equal(attr(result_double, "n_retained"), 2)
})

test_that("filter_by_singleton_ratio calculates ratio_vector attribute correctly", {
  test_table <- data.frame(
    feature_id = c("F1", "F2", "F3"),
    Sample_A = c(1, 1, 98),    # 2 singletons out of 100 = 0.02
    Sample_B = c(1, 1, 1, 97), # Oops, wrong length - fix below
    stringsAsFactors = FALSE
  )

  # Fix: proper test data
  test_table <- data.frame(
    feature_id = c("F1", "F2", "F3"),
    Sample_A = c(1, 1, 98),    # 2 singletons out of 100 = 0.02
    Sample_B = c(1, 1, 98),    # same
    stringsAsFactors = FALSE
  )

  result <- filter_by_singleton_ratio(test_table, max_singleton_ratio = 0.1)

  ratio_vec <- attr(result, "ratio_vector")
  expect_named(ratio_vec, c("Sample_A", "Sample_B"))
  expect_equal(ratio_vec["Sample_A"], 0.02, tolerance = 0.001)
  expect_equal(ratio_vec["Sample_B"], 0.02, tolerance = 0.001)
})

test_that("filter_by_singleton_ratio validates inputs", {
  test_table <- data.frame(
    feature_id = c("F1", "F2"),
    Sample_A = c(10, 20),
    stringsAsFactors = FALSE
  )

  # Invalid max_singleton_ratio
  expect_error(filter_by_singleton_ratio(test_table, max_singleton_ratio = -0.1))
  expect_error(filter_by_singleton_ratio(test_table, max_singleton_ratio = 1.5))
  expect_error(filter_by_singleton_ratio(test_table, max_singleton_ratio = "0.1"))

  # Invalid count_type
  expect_error(filter_by_singleton_ratio(test_table, count_type = "tripleton"))

  # Invalid table structure
  expect_error(filter_by_singleton_ratio(data.frame(a = 1:5)))  # No sample columns
})

test_that("filter_by_singleton_ratio handles edge cases", {
  # Table with only one sample that gets filtered
  test_table <- data.frame(
    feature_id = c("F1", "F2", "F3", "F4", "F5", "F6", "F7"),
    Bad_Sample = c(1, 1, 1, 1, 1, 1, 93),  # 6/100 = 0.06 singleton ratio
    stringsAsFactors = FALSE
  )

  result <- filter_by_singleton_ratio(test_table, max_singleton_ratio = 0.05)
  expect_equal(ncol(result), 1)  # Only feature_id column remains
  expect_equal(attr(result, "n_filtered_out"), 1)
  expect_equal(attr(result, "n_retained"), 0)
})
