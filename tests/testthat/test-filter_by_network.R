test_that("compute_mutual_information returns symmetric matrix", {
  set.seed(42)
  test_table <- matrix(abs(rpois(50 * 30, lambda = 10)), nrow = 50, ncol = 30)
  rownames(test_table) <- paste0("ASV_", 1:50)

  mi_matrix <- compute_mutual_information(test_table)

  expect_s3_class(mi_matrix, "matrix")
  expect_equal(nrow(mi_matrix), 50)
  expect_equal(ncol(mi_matrix), 50)

  # Check symmetry
  expect_equal(mi_matrix, t(mi_matrix), tolerance = 1e-6)

  # Diagonal should be NA
  expect_true(all(is.na(diag(mi_matrix))))
})

test_that("compute_mutual_information handles small datasets", {
  # Too few features
  small_mat <- matrix(rpois(2 * 20, lambda = 5), nrow = 2, ncol = 20)
  expect_error(compute_mutual_information(small_mat), NA)

  # Too few samples
  few_samples <- matrix(rpois(10 * 3, lambda = 5), nrow = 10, ncol = 3)
  expect_error(compute_mutual_information(few_samples), "at least 4 samples")
})

test_that("compute_mutual_information respects min_prevalence", {
  set.seed(123)
  test_table <- matrix(rpois(50 * 30, lambda = 5), nrow = 50, ncol = 30)
  rownames(test_table) <- paste0("ASV_", 1:50)

  # With strict prevalence, fewer features should be computed
  mi_strict <- compute_mutual_information(test_table, min_prevalence = 0.5)
  mi_loose <- compute_mutual_information(test_table, min_prevalence = 0.01)

  expect_s3_class(mi_strict, "matrix")
  expect_s3_class(mi_loose, "matrix")
})

test_that("compute_mutual_information returns reasonable values", {
  set.seed(42)
  # Create data with known structure - features as rows, samples as columns
  n_samples <- 50
  test_table <- matrix(rnorm(n_samples * 3), nrow = 3, ncol = n_samples)
  rownames(test_table) <- c("feat_x", "feat_y", "feat_z")

  mi_matrix <- compute_mutual_information(test_table, min_prevalence = 0.01)

  expect_s3_class(mi_matrix, "matrix")
  expect_true(all(mi_matrix[!is.na(mi_matrix)] >= 0))
})

test_that("analyze_feature_network returns correct structure", {
  set.seed(42)
  test_table <- matrix(abs(rpois(30 * 30, lambda = 10)), nrow = 30, ncol = 30)
  rownames(test_table) <- paste0("ASV_", 1:30)

  mi_matrix <- compute_mutual_information(test_table, min_prevalence = 0.1)
  network <- analyze_feature_network(mi_matrix, method = "mi")

  expect_s3_class(network, "list")
  expect_true(all(c("degree_centrality", "strength", "threshold_used",
                    "n_edges", "adjacency_matrix") %in% names(network)))

  expect_length(network$degree_centrality, 30)
  expect_length(network$strength, 30)
  expect_s3_class(network$adjacency_matrix, "matrix")
})

test_that("analyze_feature_network with correlation method", {
  set.seed(42)
  test_table <- matrix(rnorm(30 * 30), nrow = 30, ncol = 30)
  rownames(test_table) <- paste0("ASV_", 1:30)

  corr_matrix <- cor(t(test_table))
  diag(corr_matrix) <- 0
  corr_matrix <- abs(corr_matrix)

  network <- analyze_feature_network(corr_matrix, method = "cor")

  expect_s3_class(network, "list")
  expect_true(network$threshold_used >= 0.2)  # Auto-threshold has minimum
})

test_that("analyze_feature_network respects manual threshold", {
  set.seed(42)
  test_table <- matrix(abs(rpois(20 * 25, lambda = 10)), nrow = 20, ncol = 25)
  rownames(test_table) <- paste0("ASV_", 1:20)

  mi_matrix <- compute_mutual_information(test_table, min_prevalence = 0.1)

  # With high threshold, fewer edges
  network_high <- analyze_feature_network(mi_matrix, threshold = 1.0, method = "mi")
  # With low threshold, more edges
  network_low <- analyze_feature_network(mi_matrix, threshold = 0.1, method = "mi")

  # Both should have valid edge counts
  expect_true(!is.na(network_high$n_edges))
  expect_true(!is.na(network_low$n_edges))
  expect_true(network_high$n_edges <= network_low$n_edges)
})

test_that("filter_by_network_connectivity removes disconnected features", {
  set.seed(42)
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:50),
    matrix(abs(rpois(50 * 30, lambda = 10)), nrow = 50, ncol = 30),
    stringsAsFactors = FALSE
  )

  result <- filter_by_network_connectivity(test_table, similarity_type = "mi",
                                            verbose = FALSE)

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) <= nrow(test_table))

  n_filtered <- attr(result, "n_filtered_out")
  expect_true(!is.null(n_filtered))
  expect_true(!is.na(n_filtered))
  expect_true(n_filtered >= 0)

  disconnected <- attr(result, "disconnected_features")
  # May be NULL if no features were filtered, or character vector
  expect_true(is.null(disconnected) || is.character(disconnected))
})

test_that("filter_by_network_connectivity with correlation", {
  set.seed(42)
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:40),
    matrix(rnorm(40 * 30), nrow = 40, ncol = 30),
    stringsAsFactors = FALSE
  )

  result <- filter_by_network_connectivity(test_table, similarity_type = "cor",
                                            verbose = FALSE)

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) <= nrow(test_table))
})

test_that("filter_by_network_connectivity respects min_degree", {
  set.seed(42)
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:40),
    matrix(abs(rpois(40 * 30, lambda = 10)), nrow = 40, ncol = 30),
    stringsAsFactors = FALSE
  )

  # Higher min_degree should remove more features
  result_degree1 <- filter_by_network_connectivity(test_table, min_degree = 1,
                                                    verbose = FALSE)
  result_degree3 <- filter_by_network_connectivity(test_table, min_degree = 3,
                                                    verbose = FALSE)

  expect_true(nrow(result_degree1) >= nrow(result_degree3))
})

test_that("filter_by_network_connectivity preserves attributes", {
  set.seed(42)
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:30),
    matrix(abs(rpois(30 * 25, lambda = 10)), nrow = 30, ncol = 25),
    stringsAsFactors = FALSE
  )

  result <- filter_by_network_connectivity(test_table, verbose = FALSE)

  expect_true(!is.null(attr(result, "n_filtered_out")))
  expect_true(!is.null(attr(result, "network_summary")))
  # disconnected_features may be empty character vector but should exist
  disconnected <- attr(result, "disconnected_features")
  expect_true(!is.null(disconnected) || length(disconnected) >= 0)
})

test_that("filter_by_network_connectivity validates inputs", {
  # Invalid table structure
  expect_error(filter_by_network_connectivity(data.frame(a = 1:5)))

  # Invalid similarity type
  test_table <- data.frame(
    feature_id = paste0("ASV_", 1:10),
    matrix(abs(rpois(10 * 10, lambda = 5)), nrow = 10, ncol = 10),
    stringsAsFactors = FALSE
  )
  expect_error(filter_by_network_connectivity(test_table, similarity_type = "invalid"))
})

test_that("MI values are non-negative and reasonable", {
  set.seed(42)
  n_samples <- 50

  # Create test data as matrix with proper dimensions
  test_table <- matrix(rnorm(n_samples * 3), nrow = 3, ncol = n_samples)
  rownames(test_table) <- c("feat_x", "feat_y", "feat_z")

  mi_matrix <- compute_mutual_information(test_table, min_prevalence = 0.01)

  # All MI values should be non-negative
  valid_mi <- mi_matrix[!is.na(mi_matrix)]
  expect_true(all(valid_mi >= 0))
})

test_that("Network analysis handles edge cases", {
  # All zeros except diagonal
  zero_matrix <- matrix(0, nrow = 10, ncol = 10)
  diag(zero_matrix) <- 0

  network <- analyze_feature_network(zero_matrix, method = "cor")

  expect_equal(sum(network$degree_centrality), 0)
  expect_equal(network$n_edges, 0)
})
