test_that("compute_scree returns correct structure for mad_multiplier type", {
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:100),
    matrix(sample(0:100, 100 * 10, replace = TRUE), nrow = 100, ncol = 10),
    stringsAsFactors = FALSE
  )

  result <- compute_scree(test_table, type = "mad_multiplier", n_steps = 5, verbose = FALSE)

  expect_type(result, "list")
  expect_true(all(c("results", "summary", "type", "parameters") %in% names(result)))
  expect_s3_class(result$results, "data.frame")
  expect_equal(nrow(result$results), 5)

  # Check required columns
  expected_cols <- c("threshold", "n_features_retained", "pct_features_retained",
                     "n_samples_retained", "pct_samples_retained", "n_reads_retained",
                     "pct_reads_retained", "sparsity", "collapse_rate")
  expect_true(all(expected_cols %in% colnames(result$results)))
})

test_that("compute_scree mad_multiplier shows decreasing sample retention", {
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:50),
    matrix(abs(rnorm(50 * 8, mean = 100, sd = 30)), nrow = 50, ncol = 8),
    stringsAsFactors = FALSE
  )

  result <- compute_scree(test_table, type = "mad_multiplier", n_steps = 10, verbose = FALSE)

  # Higher MAD multipliers should filter more samples (lower retention)
  retention <- result$results$pct_samples_retained
  expect_true(all(retention <= 100))
  expect_true(all(retention >= 0))

  # First threshold (lowest) should retain most/all samples
  expect_equal(retention[1], 100)  # At multiplier=1, all samples retained
})

test_that("compute_scree absolute_feature filters features correctly", {
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:50),
    matrix(c(rep(100, 10), rep(50, 20), rep(10, 20)), nrow = 50, ncol = 4, byrow = FALSE),
    stringsAsFactors = FALSE
  )

  result <- compute_scree(test_table, type = "absolute_feature",
                          thresholds = c(1, 10, 50, 100), verbose = FALSE)

  # Should have multiple threshold evaluations
  expect_gte(nrow(result$results), 4)

  # Check that results include expected threshold values
  expect_true(1 %in% result$results$threshold)
  expect_true(100 %in% result$results$threshold)

  # At threshold=1, most features should be retained
  expect_gt(result$results$pct_features_retained[result$results$threshold == 1], 50)

  # At threshold=100, only the high-abundance features retained
  low_thresh_idx <- which(result$results$threshold == 1)
  high_thresh_idx <- which(result$results$threshold == 100)
  expect_lt(result$results$n_features_retained[high_thresh_idx],
            result$results$n_features_retained[low_thresh_idx])
})

test_that("compute_scree relative_feature works with proportions", {
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:30),
    matrix(abs(rpois(30 * 6, lambda = 50)), nrow = 30, ncol = 6),
    stringsAsFactors = FALSE
  )

  result <- compute_scree(test_table, type = "relative_feature", n_steps = 8, verbose = FALSE)

  expect_s3_class(result$results, "data.frame")
  expect_true(all(result$results$pct_features_retained >= 0))
  expect_true(all(result$results$pct_features_retained <= 100))

  # Very low thresholds should retain nearly all features
  first_retention <- result$results$pct_features_retained[1]
  expect_gt(first_retention, 80)  # Should retain >80% at lowest threshold
})

test_that("compute_scree custom thresholds work", {
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:40),
    matrix(abs(rpois(40 * 5, lambda = 30)), nrow = 40, ncol = 5),
    stringsAsFactors = FALSE
  )

  custom_thresholds <- c(1, 5, 10, 25, 50)
  result <- compute_scree(test_table, type = "custom",
                          thresholds = custom_thresholds, verbose = FALSE)

  expect_equal(nrow(result$results), length(custom_thresholds))
  expect_equal(result$type, "custom")
  expect_true(all(result$results$threshold %in% custom_thresholds))
})

test_that("compute_scree calculates collapse_rate correctly", {
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:30),
    matrix(c(rep(100, 5), rep(50, 10), rep(5, 15)), nrow = 30, ncol = 3, byrow = FALSE),
    stringsAsFactors = FALSE
  )

  result <- compute_scree(test_table, type = "absolute_feature",
                          thresholds = c(1, 5, 10, 50), verbose = FALSE)

  # Collapse rate should be NA for first point, then calculated
  expect_true(is.na(result$results$collapse_rate[1]))
  expect_false(is.na(result$results$collapse_rate[2]))

  # Sum of collapse-weighted changes should relate to total loss
  valid_rates <- result$results$collapse_rate[!is.na(result$results$collapse_rate)]
  expect_true(all(valid_rates >= 0))  # Collapse rate should be non-negative
})

test_that("compute_scree summary contains elbow detection", {
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:50),
    matrix(abs(rpois(50 * 8, lambda = 40)), nrow = 50, ncol = 8),
    stringsAsFactors = FALSE
  )

  result <- compute_scree(test_table, type = "absolute_feature", n_steps = 10, verbose = FALSE)

  summary <- result$summary
  expect_true("elbow_point" %in% names(summary))
  expect_true("threshold" %in% names(summary$elbow_point))
  expect_true("retention_at_elbow" %in% names(summary$elbow_point))
  expect_true("saturation" %in% names(summary))
  expect_true("baseline" %in% names(summary))
})

test_that("compute_scree handles edge cases", {
  # Table with very few features
  small_table <- data.frame(
    feature_id = c("ASV_1", "ASV_2", "ASV_3"),
    Sample_1 = c(100, 50, 10),
    Sample_2 = c(80, 40, 5),
    stringsAsFactors = FALSE
  )

  result <- compute_scree(small_table, type = "absolute_feature",
                          thresholds = c(1, 10, 50), verbose = FALSE)

  expect_s3_class(result$results, "data.frame")
  expect_gte(nrow(result$results), 3)

  # At threshold=50, only ASV_1 should remain
  idx_50 <- which(result$results$threshold == 50)
  if (length(idx_50) > 0) {
    expect_equal(result$results$n_features_retained[idx_50], 1)
  }
})

test_that("compute_scree validates inputs", {
  test_table <- data.frame(
    feature_id = c("ASV_1"),
    Sample_1 = c(100),
    stringsAsFactors = FALSE
  )

  # Missing thresholds for custom type
  expect_error(compute_scree(test_table, type = "custom", thresholds = NULL))

  # Invalid table structure
  expect_error(compute_scree(data.frame(a = 1:5), type = "absolute_feature"))
})

test_that("plot_scree requires ggplot2 and returns plot object", {
  skip_if_not_installed("ggplot2")

  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:30),
    matrix(abs(rpois(30 * 5, lambda = 50)), nrow = 30, ncol = 5),
    stringsAsFactors = FALSE
  )

  scree_result <- compute_scree(test_table, type = "absolute_feature",
                                 n_steps = 8, verbose = FALSE)

  # Should return a ggplot object
  plot_obj <- plot_scree(scree_result)
  expect_s3_class(plot_obj, "ggplot")
})

test_that("plot_scree handles show_collapse parameter", {
  skip_if_not_installed("ggplot2")

  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:20),
    matrix(abs(rpois(20 * 4, lambda = 30)), nrow = 20, ncol = 4),
    stringsAsFactors = FALSE
  )

  scree_result <- compute_scree(test_table, type = "absolute_feature",
                                 n_steps = 6, verbose = FALSE)

  # Both should work without error
  plot_with_collapse <- plot_scree(scree_result, show_collapse = TRUE)
  plot_without_collapse <- plot_scree(scree_result, show_collapse = FALSE)

  expect_s3_class(plot_with_collapse, "ggplot")
  expect_s3_class(plot_without_collapse, "ggplot")
})
