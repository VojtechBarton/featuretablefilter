## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 6
)

## ----installation, eval = FALSE-----------------------------------------------
# # Install from local source
# devtools::install(".")
# 
# # Required packages
# library(featuretablefilter)

## ----load-data, eval = FALSE--------------------------------------------------
# # Load TSV file (auto-detects format)
# table <- load_feature_table("feature-table.tsv")
# 
# # Load CSV file
# table <- load_feature_table("feature-table.csv")
# 
# # The first column should contain feature IDs
# head(table)

## ----coverage-absolute, eval = FALSE------------------------------------------
# # Filter samples with < 1000 reads
# filtered_table <- filter_by_coverage(table, min_reads = 1000)

## ----mad-cutoff, eval = FALSE-------------------------------------------------
# # Estimate cutoff using MAD method (default multiplier = 3)
# est <- estimate_mad_cutoff(table, multiplier = 3, floor = 0)
# print(est)
# 
# # Use the estimated cutoff for filtering
# filtered_table <- filter_by_coverage(table, min_reads = est$cutoff)

## ----iqr-cutoff, eval = FALSE-------------------------------------------------
# # Estimate cutoff using IQR method (default multiplier = 1.5)
# est <- estimate_iqr_cutoff(table, multiplier = 1.5, floor = 0)
# print(est)
# 
# # Use the estimated cutoff for filtering
# filtered_table <- filter_by_coverage(table, min_reads = est$cutoff)

## ----good-coverage, eval = FALSE----------------------------------------------
# # Estimate Good's coverage for all samples
# good_cov <- estimate_good_coverage(table, target_coverage = 0.95)
# print(good_cov)
# 
# # Filter to keep only samples with >= 95% estimated coverage
# result <- filter_by_coverage_estimator(
#   table,
#   method = "good",
#   target_coverage = 0.95
# )
# filtered_table <- result$table
# 
# # View coverage before and after filtering
# cat(sprintf("Mean coverage before: %.2f%%\n", result$mean_coverage_before * 100))
# cat(sprintf("Mean coverage after: %.2f%%\n", result$mean_coverage_after * 100))

## ----chao-coverage, eval = FALSE----------------------------------------------
# # Estimate Chao's coverage for all samples
# chao_cov <- estimate_chao_coverage(table, target_coverage = 0.90)
# print(chao_cov)
# 
# # Filter to keep only samples with >= 90% estimated coverage
# result <- filter_by_coverage_estimator(
#   table,
#   method = "chao",
#   target_coverage = 0.90
# )
# filtered_table <- result$table

## ----compare-coverage, eval = FALSE-------------------------------------------
# # Compare different estimation methods
# good_est <- estimate_good_coverage(table, target_coverage = 0.95)
# chao_est <- estimate_chao_coverage(table, target_coverage = 0.90)
# 
# cat(sprintf("Good's coverage: mean=%.2f%%, min=%.2f%%\n",
#             good_est$mean_coverage * 100, good_est$min_coverage * 100))
# cat(sprintf("Chao's coverage: mean=%.2f%%, min=%.2f%%\n",
#             chao_est$mean_coverage * 100, chao_est$min_coverage * 100))

## ----feature-absolute, eval = FALSE-------------------------------------------
# # Filter features with < 5 total reads
# filtered_table <- filter_features_by_abundance(
#   table,
#   threshold = 5,
#   mode = "absolute",
#   min_samples = 1
# )
# 
# # Keep only features present in at least 3 samples
# filtered_table <- filter_features_by_abundance(
#   table,
#   threshold = 1,
#   mode = "absolute",
#   min_samples = 3
# )

## ----feature-relative, eval = FALSE-------------------------------------------
# # Filter features with < 0.1% relative abundance
# filtered_table <- filter_features_by_abundance(
#   table,
#   threshold = 0.001,
#   mode = "relative",
#   min_samples = 1
# )

## ----relative-cutoff, eval = FALSE--------------------------------------------
# filtered_table <- filter_by_relative_cutoff(
#   table,
#   min_coverage = 1000,  # Minimum coverage required for a sample
#   relative_threshold = 0.01,  # 1% of min-coverage sample
#   remove_features = TRUE
# )

## ----pipeline, eval = FALSE---------------------------------------------------
# result <- run_filtering_pipeline(
#   input_file = "feature-table.tsv",
#   output_dir = "results",
#   prefix = "analysis1",
# 
#   # Coverage filtering
#   cov_filter_method = "mad",
#   cov_threshold = 3,
# 
#   # Abundance filtering
#   abun_filter_method = "relative_cutoff",
#   abun_threshold = 0.01,
#   min_coverage_for_relative = 1000,
# 
#   # Options
#   generate_plots = TRUE,
#   generate_report = TRUE,
#   verbose = TRUE
# )
# 
# # Access results
# original <- result$original_table
# filtered <- result$filtered_table
# qc <- result$qc_metrics

## ----qc-metrics, eval = FALSE-------------------------------------------------
# qc <- compute_filtering_qc(original_table, filtered_table, top_n = 10)
# 
# # View results
# print(qc)

## ----plot-coverage, eval = FALSE----------------------------------------------
# plot_coverage_histogram(table, threshold = 1000)

## ----plot-presence, eval = FALSE----------------------------------------------
# # Single table analysis
# plot_presence_analysis(table, threshold = 1)
# 
# # Comparison of original vs filtered
# plot_presence_analysis(
#   original_table,
#   table_filtered = filtered_table,
#   save_dir = "plots",
#   prefix = "presence_comparison"
# )

## ----plot-qc, eval = FALSE----------------------------------------------------
# plot_qc_comparison(
#   original_table,
#   filtered_table,
#   plot_dir = "plots",
#   prefix = "qc_comparison",
#   top_n = 10,
#   heatmap_top = 50
# )

## ----plot-top-features, eval = FALSE------------------------------------------
# plot_top_features_stacked(
#   original_table,
#   filtered_table,
#   top_n = 10,
#   plot_dir = "plots",
#   prefix = "top_features"
# )

## ----session-info-------------------------------------------------------------
sessionInfo()

