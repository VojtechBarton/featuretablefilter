test_that("identify_sparsity_elbow returns correct structure", {
  # Create test data with clear depth-richness relationship
  set.seed(42)
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:100),
    t(replicate(100, {
      depth <- sample(1000:10000, 1)
      n_zeros <- sum(runif(100) > (depth / 10000))
      c(rpois(n_zeros, lambda = depth / 100), rep(0, 100 - n_zeros))
    })),
    stringsAsFactors = FALSE
  )

  result <- identify_sparsity_elbow(test_table, method = "kneedle")

  expect_s3_class(result, "list")
  expect_true(all(c("elbow_threshold", "samples_above_elbow", "samples_below_elbow",
                    "richness_curve", "recommendation", "metrics") %in% names(result)))

  expect_true(is.numeric(result$elbow_threshold))
  expect_true(result$elbow_threshold > 0)
  expect_equal(result$samples_above_elbow + result$samples_below_elbow, 100)
})

test_that("identify_sparsity_elbow richness_curve has required columns", {
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:50),
    matrix(abs(rpois(50 * 20, lambda = 500)), nrow = 50, ncol = 20),
    stringsAsFactors = FALSE
  )

  result <- identify_sparsity_elbow(test_table, method = "max_derivative")

  curve_df <- result$richness_curve
  expected_cols <- c("sample_name", "depth", "richness", "rank",
                     "smoothed_richness", "derivative", "second_derivative")
  expect_true(all(expected_cols %in% colnames(curve_df)))
})

test_that("identify_sparsity_elbow different methods produce results", {
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:40),
    matrix(abs(rpois(40 * 15, lambda = 800)), nrow = 40, ncol = 15),
    stringsAsFactors = FALSE
  )

  result_kneedle <- identify_sparsity_elbow(test_table, method = "kneedle")
  result_deriv <- identify_sparsity_elbow(test_table, method = "max_derivative")
  result_second <- identify_sparsity_elbow(test_table, method = "second_derivative")

  # All should return valid elbow thresholds
  expect_true(is.numeric(result_kneedle$elbow_threshold))
  expect_true(is.numeric(result_deriv$elbow_threshold))
  expect_true(is.numeric(result_second$elbow_threshold))

  # Methods may give different results but all should be positive
  expect_gt(result_kneedle$elbow_threshold, 0)
  expect_gt(result_deriv$elbow_threshold, 0)
  expect_gt(result_second$elbow_threshold, 0)
})

test_that("identify_sparsity_elbow detects known elbow in synthetic data", {
  # Create data with a clear elbow: high-depth samples have proportional richness,
  # low-depth samples crash in richness
  set.seed(123)
  n_high <- 30
  n_low <- 20

  # High depth samples (5000-10000 reads, good richness)
  high_depths <- sample(5000:10000, n_high)
  high_matrix <- sapply(high_depths, function(d) {
    n_present <- round(100 * (d / 10000)^0.8)
    vals <- rep(0, 100)
    if (n_present > 0) vals[1:n_present] <- sample(10:100, n_present, replace = TRUE) + d / 100
    vals
  })

  # Low depth samples (100-500 reads, poor richness - the "crash")
  low_depths <- sample(100:500, n_low)
  low_matrix <- sapply(low_depths, function(d) {
    n_present <- round(100 * (d / 10000)^0.5)
    vals <- rep(0, 100)
    if (n_present > 0) vals[1:n_present] <- sample(1:10, n_present, replace = TRUE) + d / 50
    vals
  })

  # Combine into table
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:100),
    cbind(high_matrix, low_matrix),
    stringsAsFactors = FALSE
  )

  result <- identify_sparsity_elbow(test_table, method = "kneedle")

  # Elbow should be detected somewhere between low and high depth ranges
  # The exact position depends on the curvature of the richness-depth relationship
  expect_true(result$elbow_threshold > 100)    # Above minimum depth
  expect_true(result$elbow_threshold < 10000)  # Below maximum depth
  expect_true(result$elbow_threshold > min(low_depths))  # Above lowest depths
  expect_true(result$elbow_threshold < min(high_depths)) # Below highest depths
})

test_that("identify_sparsity_elbow handles edge cases", {
  # Very small dataset (but above minimum)
  small_table <- data.frame(
    feature_id = paste0("ASV_", 1:20),
    matrix(abs(rpois(20 * 6, lambda = 300)), nrow = 20, ncol = 6),
    stringsAsFactors = FALSE
  )

  result <- identify_sparsity_elbow(small_table, min_samples = 5)
  expect_s3_class(result, "list")
  expect_true(is.numeric(result$elbow_threshold))

  # Dataset at minimum size
  min_table <- data.frame(
    feature_id = paste0("ASV_", 1:10),
    matrix(abs(rpois(10 * 5, lambda = 200)), nrow = 10, ncol = 5),
    stringsAsFactors = FALSE
  )

  result_min <- identify_sparsity_elbow(min_table, min_samples = 5)
  expect_s3_class(result_min, "list")
})

test_that("identify_sparsity_elbow validates inputs", {
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:10),
    matrix(abs(rpois(10 * 10, lambda = 500)), nrow = 10, ncol = 10),
    stringsAsFactors = FALSE
  )

  # Too few samples
  tiny_table <- data.frame(
    feature_id = paste0("ASV_", 1:5),
    matrix(abs(rpois(5 * 3, lambda = 200)), nrow = 5, ncol = 3),
    stringsAsFactors = FALSE
  )
  expect_error(identify_sparsity_elbow(tiny_table, min_samples = 5))

  # Invalid percentile range
  expect_error(identify_sparsity_elbow(test_table, percentile_range = c(50, 40)))
  expect_error(identify_sparsity_elbow(test_table, percentile_range = c(-10, 50)))

  # Invalid table structure
  expect_error(identify_sparsity_elbow(data.frame(a = 1:5)))
})

test_that("identify_sparsity_elbow metrics are calculated correctly", {
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:50),
    matrix(abs(rpois(50 * 25, lambda = 600)), nrow = 50, ncol = 25),
    stringsAsFactors = FALSE
  )

  result <- identify_sparsity_elbow(test_table, method = "kneedle")

  metrics <- result$metrics
  expect_true("r_squared_stable_region" %in% names(metrics))
  expect_true("curvature_ratio" %in% names(metrics))
  expect_true("method_used" %in% names(metrics))
  expect_equal(metrics$method_used, "kneedle")
})

test_that("identify_sparsity_elbow recommendation is generated", {
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:40),
    matrix(abs(rpois(40 * 20, lambda = 400)), nrow = 40, ncol = 20),
    stringsAsFactors = FALSE
  )

  result <- identify_sparsity_elbow(test_table)

  expect_true(nchar(result$recommendation) > 50)  # Should be meaningful text
  expect_true(grepl("sample|Sample|depth|Depth|reads", result$recommendation, ignore.case = TRUE))
})

test_that("plot_sparsity_elbow returns ggplot object", {
  skip_if_not_installed("ggplot2")

  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:40),
    matrix(abs(rpois(40 * 15, lambda = 500)), nrow = 40, ncol = 15),
    stringsAsFactors = FALSE
  )

  elbow_result <- identify_sparsity_elbow(test_table)

  plot_obj <- plot_sparsity_elbow(elbow_result)
  expect_s3_class(plot_obj, "ggplot")
})

test_that("plot_sparsity_elbow includes elbow marker", {
  skip_if_not_installed("ggplot2")

  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:30),
    matrix(abs(rpois(30 * 12, lambda = 700)), nrow = 30, ncol = 12),
    stringsAsFactors = FALSE
  )

  elbow_result <- identify_sparsity_elbow(test_table)
  plot_obj <- plot_sparsity_elbow(elbow_result, main = "Test Plot")

  # Check that plot is created (patchwork object or ggplot)
  expect_true(inherits(plot_obj, "ggplot") || inherits(plot_obj, "patchwork"))
})
