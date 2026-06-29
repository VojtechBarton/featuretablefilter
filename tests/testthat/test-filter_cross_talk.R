test_that("filter_cross_talk removes low relative abundance reads", {
  # Create a feature table with clear leakage pattern:
  # Feature A: high in Sample_1 (1000), low "leakage" in others (1, 2)
  # Feature B: moderate everywhere (no leakage pattern)
  # Feature C: high in Sample_2 (500), leakage in Sample_1 (1)
  test_table <- data.frame(
    feature_id = c("ASV_A", "ASV_B", "ASV_C"),
    Sample_1 = c(1000, 100, 1),      # ASV_C is leakage here (1 << 0.001 * 500 = 0.5)
    Sample_2 = c(1, 100, 500),       # ASV_A is leakage here (1 < 0.001 * 1000 = 1)
    Sample_3 = c(2, 100, 5),         # ASV_A is leakage here (2 > 1, so kept at 0.001)
    stringsAsFactors = FALSE
  )

  # At 0.1% threshold (0.001):
  # ASV_A max = 1000, threshold = 1. So values < 1 are leakage (none, since min count is 1)
  # ASV_B max = 100, threshold = 0.1. All values >= 1, no leakage
  # ASV_C max = 500, threshold = 0.5. Values < 0.5 are leakage (none, min count is 1)

  result <- filter_cross_talk(test_table, max_rel_threshold = 0.001)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)  # All features kept
})

test_that("filter_cross_talk with stricter threshold detects leakage", {
  test_table <- data.frame(
    feature_id = c("ASV_A", "ASV_B", "ASV_C"),
    Sample_1 = c(1000, 100, 1),
    Sample_2 = c(1, 100, 500),
    Sample_3 = c(2, 100, 5),
    stringsAsFactors = FALSE
  )

  # At 0.05% threshold (0.0005):
  # ASV_A max = 1000, threshold = 0.5. Values 1 and 2 are > 0.5, kept
  # ASV_C max = 500, threshold = 0.25. Value 1 is > 0.25, kept

  # At 0.15% threshold (0.0015):
  # ASV_A max = 1000, threshold = 1.5. Value 1 < 1.5 is leakage, value 2 >= 1.5 kept
  # ASV_C max = 500, threshold = 0.75. Value 1 > 0.75, kept

  result <- filter_cross_talk(test_table, max_rel_threshold = 0.0015)

  # ASV_A in Sample_2 should be zeroed (1 < 1.5)
  expect_equal(result["ASV_A", "Sample_2"], 0)
  # ASV_A in Sample_3 should remain (2 >= 1.5)
  expect_equal(result["ASV_A", "Sample_3"], 2)
})

test_that("filter_cross_talk respects min_abs_cutoff", {
  test_table <- data.frame(
    feature_id = c("ASV_A", "ASV_B"),
    Sample_1 = c(1000, 100),
    Sample_2 = c(1, 100),
    Sample_3 = c(3, 100),
    stringsAsFactors = FALSE
  )

  # With threshold 0.001 and min_abs_cutoff = 2:
  # ASV_A max = 1000, rel_threshold = 1. Value 1 < 1 AND 1 < 2 -> leakage
  # ASV_A value 3: 3 > 1 AND 3 >= 2 -> kept

  result <- filter_cross_talk(test_table, max_rel_threshold = 0.001, min_abs_cutoff = 2)

  expect_equal(result["ASV_A", "Sample_2"], 0)  # Leakage (both conditions met)
  expect_equal(result["ASV_A", "Sample_3"], 3)  # Kept (above absolute cutoff)
})

test_that("filter_cross_talk mode='remove_feature' removes entire features", {
  test_table <- data.frame(
    feature_id = c("ASV_A", "ASV_B"),
    Sample_1 = c(1000, 100),
    Sample_2 = c(1, 100),
    Sample_3 = c(100, 100),
    stringsAsFactors = FALSE
  )

  # ASV_A has leakage (1 < 0.001 * 1000 = 1 is false, but let's use 0.002)
  # With 0.002 threshold: ASV_A threshold = 2, so 1 < 2 is leakage

  result <- filter_cross_talk(test_table, max_rel_threshold = 0.002, mode = "remove_feature")

  expect_equal(nrow(result), 1)  # Only ASV_B remains
  expect_equal(rownames(result)[1], "ASV_B")
})

test_that("filter_cross_talk mode='flag' returns original data", {
  test_table <- data.frame(
    feature_id = c("ASV_A", "ASV_B"),
    Sample_1 = c(1000, 100),
    Sample_2 = c(1, 100),
    stringsAsFactors = FALSE
  )

  result <- filter_cross_talk(test_table, max_rel_threshold = 0.002, mode = "flag")

  # Original data unchanged
  expect_equal(result["ASV_A", "Sample_2"], 1)
  expect_equal(attr(result, "n_leakage_zeros"), NA_integer_)
})

test_that("filter_cross_talk returns detailed info when requested", {
  test_table <- data.frame(
    feature_id = c("ASV_A", "ASV_B"),
    Sample_1 = c(1000, 100),
    Sample_2 = c(1, 100),
    Sample_3 = c(2, 100),
    stringsAsFactors = FALSE
  )

  result <- filter_cross_talk(test_table, max_rel_threshold = 0.0015, return_details = TRUE)

  # Check leakage matrix exists and has correct dimensions
  leakage_mat <- attr(result, "leakage_matrix")
  expect_s3_class(leakage_mat, "matrix")
  expect_equal(dim(leakage_mat), c(2, 3))

  # ASV_A in Sample_2 should be flagged (1 < 1.5)
  expect_true(leakage_mat["ASV_A", "Sample_2"])
  expect_false(leakage_mat["ASV_A", "Sample_3"])  # 2 >= 1.5

  # Check feature_max
  feature_max <- attr(result, "feature_max")
  expect_equal(feature_max["ASV_A"], 1000)
  expect_equal(feature_max["ASV_B"], 100)
})

test_that("filter_cross_talk handles edge cases", {
  # Table with all zeros for one feature
  test_table <- data.frame(
    feature_id = c("ASV_A", "ASV_B"),
    Sample_1 = c(1000, 0),
    Sample_2 = c(1, 0),
    stringsAsFactors = FALSE
  )

  # Should not error on zero-max features
  result <- filter_cross_talk(test_table, max_rel_threshold = 0.001)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
})

test_that("filter_cross_talk validates inputs", {
  test_table <- data.frame(
    feature_id = c("ASV_A"),
    Sample_1 = c(100),
    stringsAsFactors = FALSE
  )

  expect_error(filter_cross_talk(test_table, max_rel_threshold = -0.1))
  expect_error(filter_cross_talk(test_table, max_rel_threshold = 1.5))
  expect_error(filter_cross_talk(test_table, min_abs_cutoff = -1))
  expect_error(filter_cross_talk(data.frame(a = 1:5)))  # No sample columns
})

test_that("filter_index_hopping alias works identically", {
  test_table <- data.frame(
    feature_id = c("ASV_A", "ASV_B"),
    Sample_1 = c(1000, 100),
    Sample_2 = c(1, 100),
    stringsAsFactors = FALSE
  )

  result1 <- filter_cross_talk(test_table, max_rel_threshold = 0.002)
  result2 <- filter_index_hopping(test_table, max_rel_threshold = 0.002)

  expect_equal(result1, result2)
})
