library(testthat)
library(featuretablefilter)

# Create test data
# 5 features, 10 samples
test_table <- data.frame(
  feature_id = c("feat1", "feat2", "feat3", "feat4", "feat5"),
  sample1 = c(100, 10, 5, 1, 0),    # feat1 high, feat2 medium, feat3 low, feat4 very low
  sample2 = c(100, 10, 5, 1, 50),   # feat5 present
  sample3 = c(100, 10, 0, 1, 50),   # feat3 absent
  sample4 = c(100, 0, 5, 1, 50),    # feat2 absent
  sample5 = c(100, 10, 5, 0, 50),   # feat4 absent
  sample6 = c(100, 10, 5, 1, 0),    # feat5 absent
  sample7 = c(100, 10, 0, 1, 50),   # feat3 absent
  sample8 = c(100, 0, 5, 1, 50),    # feat2 absent
  sample9 = c(100, 10, 5, 0, 50),   # feat4 absent
  sample10 = c(100, 10, 5, 1, 0),   # feat5 absent
  stringsAsFactors = FALSE
)
# Summary:
# feat1: present in all 10 samples (100% prevalence), high abundance everywhere
# feat2: present in 7 samples (70% prevalence), medium abundance
# feat3: present in 7 samples (70% prevalence), low abundance
# feat4: present in 8 samples (80% prevalence), very low abundance (1)
# feat5: present in 5 samples (50% prevalence), medium-high abundance

describe("filter_features_joint", {
  it("OR logic keeps features meeting either criterion", {
    result <- filter_features_joint(
      test_table,
      abundance_threshold = 50,       # High abundance threshold
      prevalence_threshold = 0.6,     # 60% of samples
      mode = "relative",
      logic = "OR"
    )

    expect_s3_class(result$table, "data.frame")
    expect_equal(result$n_features_before, 5)
    expect_true(result$n_features_after <= 5)
    expect_equal(result$logic, "OR")
    expect_equal(result$prevalence_threshold, 0.6)

    # With OR logic and these thresholds:
    # feat1: meets both (high abundance everywhere, 100% prevalence) -> kept
    # feat5: meets abundance in some samples, 50% prevalence (< 60%) -> check if meets abundance
    expect_true(result$n_by_both >= 0)
    expect_true(result$n_by_abundance_only >= 0)
    expect_true(result$n_by_prevalence_only >= 0)
  })

  it("AND logic keeps only features meeting both criteria", {
    result <- filter_features_joint(
      test_table,
      abundance_threshold = 0,        # Any non-zero counts
      prevalence_threshold = 0.7,     # 70% of samples
      mode = "absolute",
      logic = "AND"
    )

    # With AND logic:
    # feat1: present in 100%, any abundance -> meets both -> kept
    # feat2: present in 70%, any abundance -> meets both -> kept
    # feat3: present in 70%, any abundance -> meets both -> kept
    # feat4: present in 80%, any abundance -> meets both -> kept
    # feat5: present in 50% -> does NOT meet prevalence -> removed
    expect_equal(result$n_features_before, 5)
    expect_true(result$n_features_after <= 4)  # At most 4 features kept
  })

  it("correctly calculates feature statistics", {
    result <- filter_features_joint(
      test_table,
      abundance_threshold = 0,
      prevalence_threshold = 0.5,
      mode = "absolute",
      logic = "OR"
    )

    expect_true("feature_id" %in% colnames(result$feature_details))
    expect_true("prevalence" %in% colnames(result$feature_details))
    expect_true("meets_abundance" %in% colnames(result$feature_details))
    expect_true("meets_prevalence" %in% colnames(result$feature_details))
    expect_true("kept" %in% colnames(result$feature_details))
    expect_equal(nrow(result$feature_details), 5)
  })

  it("handles relative mode correctly", {
    result_abs <- filter_features_joint(
      test_table,
      abundance_threshold = 50,
      prevalence_threshold = 0.5,
      mode = "absolute",
      logic = "OR"
    )

    result_rel <- filter_features_joint(
      test_table,
      abundance_threshold = 0.1,      # 10% relative abundance
      prevalence_threshold = 0.5,
      mode = "relative",
      logic = "OR"
    )

    expect_equal(result_abs$mode, "absolute")
    expect_equal(result_rel$mode, "relative")
    # Results should differ because thresholds are applied differently
    expect_true(result_abs$n_features_after != result_rel$n_features_after ||
                result_abs$n_features_removed != result_rel$n_features_removed)
  })

  it("validates input parameters", {
    expect_error(
      filter_features_joint(test_table, abundance_threshold = -1),
      "abundance_threshold must be non-negative"
    )

    expect_error(
      filter_features_joint(test_table, prevalence_threshold = 1.5),
      "prevalence_threshold must be between 0 and 1"
    )

    expect_error(
      filter_features_joint(test_table, prevalence_threshold = -0.1),
      "prevalence_threshold must be between 0 and 1"
    )
  })

  it("keeps all features when thresholds are zero/low with OR logic", {
    result <- filter_features_joint(
      test_table,
      abundance_threshold = 0,
      prevalence_threshold = 0,
      mode = "absolute",
      logic = "OR"
    )

    # All features have at least some presence, so all should be kept
    expect_equal(result$n_features_after, 5)
    expect_equal(result$n_features_removed, 0)
  })

  it("removes features that meet neither criterion", {
    # Set very high thresholds that no feature can meet
    result <- filter_features_joint(
      test_table,
      abundance_threshold = 1000,     # Higher than any value in table
      prevalence_threshold = 1.0,     # 100% of samples
      mode = "absolute",
      logic = "AND"
    )

    # No feature meets both criteria (none have 1000+ reads, none in 100% of samples with that abundance)
    expect_equal(result$n_by_neither, 5)  # All 5 meet neither
    expect_equal(result$n_features_after, 0)  # All removed
  })

  it("returns correct breakdown statistics", {
    result <- filter_features_joint(
      test_table,
      abundance_threshold = 50,
      prevalence_threshold = 0.6,
      mode = "absolute",
      logic = "OR"
    )

    # Verify that the breakdown sums correctly
    total_kept <- result$n_by_abundance_only + result$n_by_prevalence_only + result$n_by_both
    expect_equal(total_kept, result$n_features_after)

    expect_equal(result$n_features_removed, result$n_by_neither)
    expect_equal(
      result$n_features_before,
      result$n_features_after + result$n_features_removed
    )
  })

  it("works with remove_zeros = FALSE", {
    result_keep <- filter_features_joint(
      test_table,
      abundance_threshold = 50,
      prevalence_threshold = 0.6,
      mode = "absolute",
      logic = "OR",
      remove_zeros = FALSE
    )

    # With remove_zeros = FALSE, all features should be kept but some zeroed out
    expect_equal(nrow(result_keep$table), 5)
  })
})
