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

test_that("compute_filtering_qc calculates Shannon ENS", {
  # Create table with known diversity
  original <- data.frame(
    id = c("f1", "f2", "f3", "f4"),
    s1 = c(25, 25, 25, 25),  # Even distribution: high diversity
    s2 = c(25, 25, 25, 25)
  )

  # Filter removes one feature
  filtered <- data.frame(
    id = c("f1", "f2", "f3"),
    s1 = c(33.3, 33.3, 33.4),  # Still relatively even
    s2 = c(33.3, 33.3, 33.4)
  )

  qc <- compute_filtering_qc(original, filtered)

  expect_type(qc$shannon_ens_original, "double")
  expect_type(qc$shannon_ens_filtered, "double")
  expect_type(qc$shannon_ens_retention_percent, "double")

  # For perfect even distribution of 4 species, Shannon ENS = 4
  # For perfect even distribution of 3 species, Shannon ENS = 3
  expect_true(qc$shannon_ens_original > 3.5)  # Close to 4
  expect_true(qc$shannon_ens_filtered > 2.5)  # Close to 3
})

test_that("compute_filtering_qc calculates Simpson ENS", {
  # Create table with known Simpson diversity
  original <- data.frame(
    id = c("f1", "f2", "f3", "f4", "f5"),
    s1 = c(20, 20, 20, 20, 20),  # Perfectly even: Simpson ENS = 5
    s2 = c(20, 20, 20, 20, 20)
  )

  filtered <- original

  qc <- compute_filtering_qc(original, filtered)

  expect_type(qc$simpson_ens_original, "double")
  expect_type(qc$simpson_ens_filtered, "double")
  expect_type(qc$simpson_ens_retention_percent, "double")

  # For perfectly even distribution of 5 species, Simpson ENS = 5
  expect_equal(qc$simpson_ens_original, 5, tolerance = 0.1)
  expect_equal(qc$simpson_ens_retention_percent, 100, tolerance = 0.1)
})

test_that("ENS retention decreases when filtering dominant features", {
  # Original: one dominant feature, many rare
  original <- data.frame(
    id = c("dominant", "rare1", "rare2", "rare3", "rare4"),
    s1 = c(80, 5, 5, 5, 5),  # Low evenness
    s2 = c(80, 5, 5, 5, 5)
  )

  # Filter removes the dominant feature
  filtered <- data.frame(
    id = c("rare1", "rare2", "rare3", "rare4"),
    s1 = c(25, 25, 25, 25),  # Now perfectly even
    s2 = c(25, 25, 25, 25)
  )

  qc <- compute_filtering_qc(original, filtered)

  # Filtering the dominant feature should actually INCREASE diversity
  expect_true(qc$shannon_ens_retention_percent > 100)
  expect_true(qc$simpson_ens_retention_percent > 100)
})

test_that("ENS metrics handle low-diversity samples", {
  # Sample with only one feature (zero diversity)
  original <- data.frame(
    id = c("f1", "f2"),
    s1 = c(100, 0),  # Only f1 present
    s2 = c(50, 50)   # Both present
  )

  filtered <- data.frame(
    id = c("f1"),
    s1 = c(100),
    s2 = c(50)
  )

  qc <- compute_filtering_qc(original, filtered)

  # Should not error, may return NA for some metrics
  expect_type(qc$shannon_ens_original, "double")
  expect_type(qc$simpson_ens_original, "double")
})

test_that("calc_shannon_ens calculates correct values", {
  # Perfectly even distribution of 4 species should give ENS ≈ 4
  mat_even <- matrix(c(25, 25, 25, 25), nrow = 4, ncol = 1)
  ens <- calc_shannon_ens(mat_even)
  expect_equal(ens, 4, tolerance = 0.01)

  # Uneven distribution should give lower ENS
  mat_uneven <- matrix(c(70, 10, 10, 10), nrow = 4, ncol = 1)
  ens_uneven <- calc_shannon_ens(mat_uneven)
  expect_true(ens_uneven < 4)
  expect_true(ens_uneven > 1)
})

test_that("calc_simpson_ens calculates correct values", {
  # Perfectly even distribution of 5 species should give ENS = 5
  mat_even <- matrix(c(20, 20, 20, 20, 20), nrow = 5, ncol = 1)
  ens <- calc_simpson_ens(mat_even)
  expect_equal(ens, 5, tolerance = 0.01)

  # Single species should give ENS = 1
  mat_single <- matrix(c(100), nrow = 1, ncol = 1)
  ens_single <- calc_simpson_ens(mat_single)
  expect_equal(as.numeric(ens_single), 1, tolerance = 0.01)
})

test_that("calc_hill_numbers returns profile across q values", {
  mat <- matrix(c(25, 25, 25, 25), nrow = 4, ncol = 1)

  hill <- calc_hill_numbers(mat, q = c(0, 1, 2))

  expect_s3_class(hill, "matrix")
  expect_equal(ncol(hill), 3)
  expect_equal(unname(hill[, 1]), 4)  # q=0 is richness
  expect_equal(unname(hill[, 2]), 4, tolerance = 0.01)  # q=1 is Shannon ENS
  expect_equal(unname(hill[, 3]), 4, tolerance = 0.01)  # q=2 is Simpson ENS

  # For perfectly even distribution, all Hill numbers should equal S
})

test_that("calc_hill_numbers shows diversity profile for uneven community", {
  # Community dominated by one species
  mat <- matrix(c(80, 5, 5, 5, 5), nrow = 5, ncol = 1)

  hill <- calc_hill_numbers(mat, q = c(0, 1, 2))

  # q=0 (richness) counts all species
  expect_equal(unname(hill[, 1]), 5)

  # Higher q gives lower ENS as it downweights rare species
  expect_true(unname(hill[, 1]) >= unname(hill[, 2]))  # Richness >= Shannon ENS
  expect_true(unname(hill[, 2]) >= unname(hill[, 3]))  # Shannon ENS >= Simpson ENS
})
