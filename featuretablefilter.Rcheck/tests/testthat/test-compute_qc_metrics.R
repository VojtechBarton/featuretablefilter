test_that("compute_filtering_qc calculates sparsity correctly", {
  # Create simple test tables
  original <- data.frame(
    id = c("f1", "f2"),
    s1 = c(10, 0),
    s2 = c(0, 5),
    s3 = c(0, 0)
  )

  filtered <- data.frame(
    id = c("f1", "f2"),
    s1 = c(10, 0),
    s2 = c(0, 5)
  )

  qc <- compute_filtering_qc(original, filtered)

  expect_type(qc$sparsity_original, "double")
  expect_type(qc$sparsity_filtered, "double")
  expect_true(qc$sparsity_original >= 0 && qc$sparsity_original <= 1)
})

test_that("compute_filtering_qc calculates retention rates", {
  original <- data.frame(
    id = c("f1", "f2", "f3"),
    s1 = c(10, 5, 3),
    s2 = c(8, 4, 2)
  )

  filtered <- data.frame(
    id = c("f1", "f2"),
    s1 = c(10, 5),
    s2 = c(8, 4)
  )

  qc <- compute_filtering_qc(original, filtered)

  expect_type(qc$read_retention_percent, "double")
  expect_type(qc$feature_retention_percent, "double")
  expect_equal(qc$feature_retention_percent, 2/3 * 100, tolerance = 0.01)
})

test_that("compute_filtering_qc returns top N features", {
  original <- data.frame(
    id = c("f1", "f2", "f3", "f4", "f5"),
    s1 = c(100, 50, 30, 20, 10),
    s2 = c(100, 50, 30, 20, 10)
  )

  filtered <- original

  qc <- compute_filtering_qc(original, filtered, top_n = 3)

  expect_length(qc$orig_top_features, 3)
  expect_length(qc$filt_top_features, 3)
  expect_true("f1" %in% qc$orig_top_features)  # Most abundant
  expect_true("f2" %in% qc$orig_top_features)  # Second most abundant
})

test_that("compute_filtering_qc calculates overlap metrics", {
  original <- data.frame(
    id = c("f1", "f2", "f3"),
    s1 = c(100, 50, 30),
    s2 = c(100, 50, 30)
  )

  filtered <- original

  qc <- compute_filtering_qc(original, filtered, top_n = 2)

  expect_type(qc$top_n_overlap_count, "integer")
  expect_type(qc$top_n_overlap_percent, "double")
  expect_type(qc$top_n_jaccard_similarity, "double")
  expect_equal(qc$top_n_jaccard_similarity, 1.0)  # Identical tables
})
