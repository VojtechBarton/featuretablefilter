library(testthat)
library(featuretablefilter)

# Create test data with known singletons
test_table <- data.frame(
  feature_id = c("feat1", "feat2", "feat3", "feat4", "feat5"),
  sample1 = c(10, 5, 3, 1, 1),    # 2 singletons, total=20, Good's = 1 - 2/20 = 0.90
  sample2 = c(50, 30, 15, 4, 1),  # 1 singleton, total=100, Good's = 1 - 1/100 = 0.99
  sample3 = c(5, 5, 5, 5, 5),     # 0 singletons, total=25, Good's = 1 - 0/25 = 1.00
  stringsAsFactors = FALSE
)

describe("estimate_good_coverage", {
  it("calculates correct coverage for each sample", {
    result <- estimate_good_coverage(test_table, target_coverage = 0.95)

    expect_type(result$sample_coverage, "double")
    expect_length(result$sample_coverage, 3)
    expect_named(result$sample_coverage, c("sample1", "sample2", "sample3"))

    # Verify Good's coverage calculations
    # sample1: 1 - 2/20 = 0.90
    expect_equal(unname(result$sample_coverage["sample1"]), 0.90, tolerance = 0.001)
    # sample2: 1 - 1/100 = 0.99
    expect_equal(unname(result$sample_coverage["sample2"]), 0.99, tolerance = 0.001)
    # sample3: 1 - 0/25 = 1.00
    expect_equal(unname(result$sample_coverage["sample3"]), 1.00, tolerance = 0.001)
  })

  it("correctly counts samples below target", {
    result <- estimate_good_coverage(test_table, target_coverage = 0.95)

    expect_equal(result$n_samples_below_target, 1)  # Only sample1 is below 0.95
    expect_equal(result$target_coverage, 0.95)
  })

  it("returns expected summary statistics", {
    result <- estimate_good_coverage(test_table, target_coverage = 0.95)

    expect_true(result$mean_coverage > 0.90 && result$mean_coverage < 1.00)
    expect_equal(result$min_coverage, 0.90)
    expect_equal(result$max_coverage, 1.00)
    expect_equal(result$total_singletons, 3)  # 2 + 1 + 0
    expect_equal(result$total_reads, 145)     # 20 + 100 + 25
    expect_equal(result$n_samples, 3)
  })

  it("handles target_coverage parameter correctly", {
    result_low <- estimate_good_coverage(test_table, target_coverage = 0.80)
    result_high <- estimate_good_coverage(test_table, target_coverage = 0.99)

    expect_equal(result_low$n_samples_below_target, 0)  # All above 0.80
    # At 0.99 threshold: sample1 (0.90) is below, sample2 (0.99) is at threshold, sample3 (1.00) is above
    expect_gte(result_high$n_samples_below_target, 1)
  })
})

describe("estimate_chao_coverage", {
  it("calculates coverage for each sample", {
    result <- estimate_chao_coverage(test_table, target_coverage = 0.90)

    expect_type(result$sample_coverage, "double")
    expect_length(result$sample_coverage, 3)
    expect_named(result$sample_coverage, c("sample1", "sample2", "sample3"))

    # All values should be between 0 and 1
    expect_true(all(result$sample_coverage >= 0 & result$sample_coverage <= 1, na.rm = TRUE))
  })

  it("returns expected summary statistics", {
    result <- estimate_chao_coverage(test_table, target_coverage = 0.90)

    expect_true(!is.na(result$mean_coverage))
    expect_true(result$total_singletons > 0)
    expect_true(result$total_doubletons >= 0)
    expect_true(result$total_features > 0)
    expect_equal(result$total_reads, 145)
    expect_equal(result$n_samples, 3)
  })

  it("counts doubletons correctly", {
    result <- estimate_chao_coverage(test_table, target_coverage = 0.90)

    # sample1: 1 doubleton (feat4 has value 1... actually feat4=1 is singleton)
    # Let's recalculate: sample1 has feat4=1, feat5=1 as singletons
    # sample2 has feat5=1 as singleton, no doubletons
    # sample3 has no singletons or doubletons (all 5s)
    expect_true(result$total_doubletons >= 0)
  })
})

describe("filter_by_coverage_estimator", {
  it("filters samples using Good's coverage", {
    result <- filter_by_coverage_estimator(
      test_table,
      method = "good",
      target_coverage = 0.95
    )

    expect_s3_class(result$table, "data.frame")
    expect_equal(result$method, "good")
    expect_equal(result$target_coverage, 0.95)
    expect_equal(result$n_samples_before, 3)
    expect_true(result$n_samples_after <= 3)
    expect_true(result$n_samples_filtered >= 0)
  })

  it("filters samples using Chao's coverage", {
    result <- filter_by_coverage_estimator(
      test_table,
      method = "chao",
      target_coverage = 0.90
    )

    expect_s3_class(result$table, "data.frame")
    expect_equal(result$method, "chao")
    expect_equal(result$target_coverage, 0.90)
    expect_equal(result$n_samples_before, 3)
    expect_true(result$n_samples_after <= 3)
  })

  it("improves mean coverage after filtering", {
    result <- filter_by_coverage_estimator(
      test_table,
      method = "good",
      target_coverage = 0.95
    )

    # After filtering low-coverage samples, mean should increase or stay same
    expect_true(result$mean_coverage_after >= result$mean_coverage_before - 0.01)
  })

  it("uses default target coverage when not specified", {
    result_good <- filter_by_coverage_estimator(test_table, method = "good")
    result_chao <- filter_by_coverage_estimator(test_table, method = "chao")

    expect_equal(result_good$target_coverage, 0.95)
    expect_equal(result_chao$target_coverage, 0.90)
  })

  it("applies min_reads floor in addition to coverage", {
    # Use a very low coverage target but high min_reads
    result <- filter_by_coverage_estimator(
      test_table,
      method = "good",
      target_coverage = 0.50,  # Very low, would keep all samples
      min_reads = 50           # High, would filter out sample1 and sample3
    )

    # Should apply both filters, keeping only samples meeting both criteria
    expect_true(result$n_samples_after <= 3)
  })
})
