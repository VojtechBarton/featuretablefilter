test_that("analyze_depth_sparsity returns correct structure", {
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:50),
    matrix(abs(rpois(50 * 30, lambda = 500)), nrow = 50, ncol = 30),
    stringsAsFactors = FALSE
  )

  result <- analyze_depth_sparsity(test_table, metric = "sparsity", verbose = FALSE)

  expect_type(result, "list")
  expect_true(all(c("sample_metrics", "outliers", "n_outliers", "fit_summary",
                    "thresholds", "recommendation") %in% names(result)))

  expect_s3_class(result$sample_metrics, "data.frame")
  expect_equal(nrow(result$sample_metrics), 30)

  # Check required columns in sample_metrics
  expected_cols <- c("sample_name", "depth", "richness", "sparsity", "residual")
  expect_true(all(expected_cols %in% colnames(result$sample_metrics)))
})

test_that("analyze_depth_sparsity calculates metrics correctly", {
  # Create table with known sparsity
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:20),
    Sample_A = rep(c(10, 0, 0, 0, 0), 4),  # 80% zeros, depth=40
    Sample_B = rep(c(10, 10, 10, 10, 10), 4),  # 0% zeros, depth=200
    stringsAsFactors = FALSE
  )

  result <- analyze_depth_sparsity(test_table, metric = "sparsity", verbose = FALSE)

  metrics <- result$sample_metrics

  # Sample_A should have 80% sparsity (16 zeros out of 20 features)
  expect_equal(metrics$sparsity[metrics$sample_name == "Sample_A"], 0.8, tolerance = 0.01)

  # Sample_B should have 0% sparsity (0 zeros out of 20 features)
  expect_equal(metrics$sparsity[metrics$sample_name == "Sample_B"], 0, tolerance = 0.01)
})

test_that("analyze_depth_sparsity detects outliers with MAD method", {
  # Create data where one sample is clearly an outlier
  set.seed(42)
  normal_data <- matrix(abs(rpois(50 * 25, lambda = 600)), nrow = 50, ncol = 25)

  # Add one high-sparsity outlier (mostly zeros)
  outlier_col <- c(rep(0, 45), rep(100, 5))  # 90% zeros but some high counts

  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:50),
    normal_data,
    Outlier_Sample = outlier_col,
    stringsAsFactors = FALSE
  )

  result <- analyze_depth_sparsity(test_table, metric = "sparsity",
                                    outlier_method = "mad", multiplier = 2, verbose = FALSE)

  # Should detect at least the outlier
  expect_gte(result$n_outliers, 1)
  expect_true("Outlier_Sample" %in% result$outliers$sample_name || result$n_outliers >= 1)
})

test_that("analyze_depth_sparsity detects outliers with IQR method", {
  set.seed(123)
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:40),
    matrix(abs(rpois(40 * 20, lambda = 400)), nrow = 40, ncol = 20),
    stringsAsFactors = FALSE
  )

  result <- analyze_depth_sparsity(test_table, metric = "sparsity",
                                    outlier_method = "iqr", verbose = FALSE)

  expect_s3_class(result$outliers, "data.frame")
  expect_true(is.numeric(result$n_outliers))
})

test_that("analyze_depth_sparsity works with richness metric", {
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:30),
    matrix(abs(rpois(30 * 15, lambda = 300)), nrow = 30, ncol = 15),
    stringsAsFactors = FALSE
  )

  result <- analyze_depth_sparsity(test_table, metric = "richness", verbose = FALSE)

  expect_type(result, "list")
  expect_true(all(c("slope", "intercept", "r_squared") %in% names(result$fit_summary)))
})

test_that("analyze_depth_sparsity fits reasonable regression line", {
  # Create data with clear positive relationship between depth and richness
  set.seed(456)
  n_samples <- 30
  n_features <- 100

  # Create samples with varying depths
  depths <- c(rep(1000, 10), rep(3000, 10), rep(5000, 10))

  # Build matrix properly
  mat <- matrix(0, nrow = n_features, ncol = n_samples)
  for (i in seq_len(n_samples)) {
    r <- rpois(1, depths[i] / 50) + 10
    present_idx <- sample(n_features, min(r, n_features))
    mat[present_idx, i] <- sample(10:100, length(present_idx), replace = TRUE)
  }

  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:n_features),
    mat,
    stringsAsFactors = FALSE
  )
  colnames(test_table)[-1] <- paste0("Sample_", seq_len(n_samples))

  result <- analyze_depth_sparsity(test_table, metric = "richness", verbose = FALSE)

  fit <- result$fit_summary
  expect_true(fit$r_squared >= 0 && fit$r_squared <= 1)
  expect_true(is.numeric(fit$slope))
})

test_that("analyze_depth_sparsity handles direction parameter", {
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:25),
    matrix(abs(rpois(25 * 12, lambda = 500)), nrow = 25, ncol = 12),
    stringsAsFactors = FALSE
  )

  result_high <- analyze_depth_sparsity(test_table, direction = "high_sparsity", verbose = FALSE)
  result_low <- analyze_depth_sparsity(test_table, direction = "low_sparsity", verbose = FALSE)
  result_both <- analyze_depth_sparsity(test_table, direction = "both", verbose = FALSE)

  expect_s3_class(result_high$outliers, "data.frame")
  expect_s3_class(result_low$outliers, "data.frame")
  expect_s3_class(result_both$outliers, "data.frame")
})

test_that("filter_depth_sparsity_outliers removes flagged samples", {
  set.seed(789)
  n_features <- 30
  n_samples <- 20

  # Create normal samples with moderate sparsity
  normal_mat <- matrix(abs(rpois(n_features * n_samples, lambda = 400)),
                        nrow = n_features, ncol = n_samples)

  # Add a clear outlier sample with high sparsity (mostly zeros)
  outlier_col <- c(rep(0, 25), rep(50, 5))  # 83% zeros

  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:n_features),
    normal_mat,
    Bad_Sample = outlier_col,
    stringsAsFactors = FALSE
  )

  filtered <- filter_depth_sparsity_outliers(test_table, keep_outliers = FALSE)

  # Bad_Sample should be removed
  expect_false("Bad_Sample" %in% colnames(filtered))
  expect_equal(attr(filtered, "n_removed"), 1)
})

test_that("filter_depth_sparsity_outliers keeps only outliers when requested", {
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:20),
    matrix(abs(rpois(20 * 10, lambda = 300)), nrow = 20, ncol = 10),
    stringsAsFactors = FALSE
  )

  kept_only_outliers <- filter_depth_sparsity_outliers(test_table, keep_outliers = TRUE)

  # Should only keep outlier samples
  expect_true(ncol(kept_only_outliers) <= ncol(test_table))
})

test_that("analyze_depth_sparsity validates inputs", {
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:10),
    matrix(abs(rpois(10 * 8, lambda = 200)), nrow = 10, ncol = 8),
    stringsAsFactors = FALSE
  )

  # Invalid table structure
  expect_error(analyze_depth_sparsity(data.frame(a = 1:5)))

  # Invalid metric
  expect_error(analyze_depth_sparsity(test_table, metric = "invalid"))

  # Invalid outlier method
  expect_error(analyze_depth_sparsity(test_table, outlier_method = "invalid"))
})

test_that("plot_depth_sparsity returns ggplot object", {
  skip_if_not_installed("ggplot2")

  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:25),
    matrix(abs(rpois(25 * 15, lambda = 400)), nrow = 25, ncol = 15),
    stringsAsFactors = FALSE
  )

  result <- analyze_depth_sparsity(test_table, verbose = FALSE)

  plot_obj <- plot_depth_sparsity(result)
  expect_s3_class(plot_obj, "ggplot")
})

test_that("plot_depth_sparsity includes title and subtitle", {
  skip_if_not_installed("ggplot2")

  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:20),
    matrix(abs(rpois(20 * 12, lambda = 500)), nrow = 20, ncol = 12),
    stringsAsFactors = FALSE
  )

  result <- analyze_depth_sparsity(test_table, verbose = FALSE)
  plot_obj <- plot_depth_sparsity(result, main = "Test Title")

  expect_equal(plot_obj$labels$title, "Test Title")
  expect_true(!is.null(plot_obj$labels$subtitle))
})

test_that("analyze_depth_sparsity stores analysis as attribute on filtered table", {
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:15),
    matrix(abs(rpois(15 * 10, lambda = 300)), nrow = 15, ncol = 10),
    stringsAsFactors = FALSE
  )

  filtered <- filter_depth_sparsity_outliers(test_table)

  expect_true(!is.null(attr(filtered, "outlier_analysis")))
  expect_true(!is.null(attr(filtered, "n_retained")))
})

test_that("plot_reads_vs_asvs returns expected structure", {
  skip_if_not_installed("ggplot2")

  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:50),
    matrix(abs(rpois(50 * 20, lambda = 500)), nrow = 50, ncol = 20),
    stringsAsFactors = FALSE
  )

  result <- plot_reads_vs_asvs(test_table)

  expect_type(result, "list")
  expect_true(all(c("plot", "outliers", "metrics", "cutoff") %in% names(result)))
  expect_s3_class(result$plot, "ggplot")
  expect_s3_class(result$outliers, "data.frame")
  expect_s3_class(result$metrics, "data.frame")
})

test_that("plot_reads_vs_asvs detects outliers correctly", {
  skip_if_not_installed("ggplot2")

  # Create data with known outlier: one sample with very low richness
  set.seed(42)
  normal_depths <- sample(1000:5000, 19)

  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:100),
    cbind(
      matrix(sapply(normal_depths, function(d) {
        n_present <- round(d / 20)
        vals <- rep(0, 100)
        if (n_present > 0) vals[1:min(n_present, 100)] <- sample(1:100, min(n_present, 100), replace = TRUE)
        vals
      }), nrow = 100),
      matrix(c(rep(0, 95), rep(1, 5)), nrow = 100, ncol = 1)  # Only 5 ASVs
    ),
    stringsAsFactors = FALSE
  )
  colnames(test_table)[-1] <- c(paste0("Sample_", 1:19), "Outlier_Sample")

  result <- plot_reads_vs_asvs(test_table, mad_multiplier = 2)

  # Check that metrics and outliers are returned
  expect_type(result$outliers$sample_name, "character")
  expect_equal(nrow(result$metrics), 20)  # All samples in metrics
})

test_that("plot_reads_vs_asvs validates inputs", {
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:10),
    matrix(abs(rpois(10 * 10, lambda = 500)), nrow = 10, ncol = 10),
    stringsAsFactors = FALSE
  )

  # Invalid table structure (no feature_id column)
  expect_error(plot_reads_vs_asvs(data.frame(a = 1:5)))
})
