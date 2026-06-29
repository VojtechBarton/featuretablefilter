#' Run a complete feature table filtering pipeline with QC report
#'
#' A comprehensive pipeline that loads, filters, and analyzes a feature table.
#' Allows selection of different coverage and abundance filtering methods.
#' Produces QC metrics and visualizations comparing before/after filtering.
#'
#' @param input_file Path to input feature table file (TSV or CSV).
#' @param output_dir Directory to save outputs (filtered table, plots, report). Default is current directory.
#' @param prefix Prefix for output filenames. Default is "filtered".
#'
#' # Coverage Filtering Options
#' @param cov_filter_method Method for sample coverage filtering:
#'                          "none" - no coverage filtering
#'                          "absolute" - use fixed min_reads threshold
#'                          "mad" - use MAD-based threshold estimation
#'                          "iqr" - use IQR/Tukey-based threshold estimation
#'                          "good" - use Good's coverage estimator (ecological completeness)
#'                          "chao" - use Chao's coverage estimator (conservative completeness)
#' @param cov_threshold Fixed threshold if method="absolute", or multiplier if method="mad"/"iqr".
#'                      For MAD/IQR, this multiplies the spread statistic.
#'                      Not used for "good" or "chao" methods (use cov_target_coverage instead).
#' @param cov_floor Minimum possible coverage cutoff when using MAD/IQR methods.
#' @param cov_target_coverage Target ecological coverage for "good" or "chao" methods (0-1).
#'                            Default is 0.95 for Good's, 0.90 for Chao's.
#' @param cov_min_reads Optional minimum absolute read count applied in addition to coverage estimators.
#'                      Useful as a safety floor. Default is 0 (no additional cutoff).
#'
#' # Abundance Filtering Options
#' @param abun_filter_method Method for feature abundance filtering:
#'                           "none" - no abundance filtering
#'                           "absolute" - filter features below absolute count threshold
#'                           "relative" - filter features below relative abundance threshold
#'                           "relative_cutoff" - use relative threshold based on min-coverage sample
#' @param abun_threshold Threshold value for abundance filtering.
#'                       For "absolute": minimum read count per feature.
#'                       For "relative": proportion (e.g., 0.01 for 1%).
#'                       For "relative_cutoff": relative proportion used to calculate absolute threshold.
#' @param abun_min_samples Minimum number of samples where feature must exceed threshold (default: 1).
#'
#' # Additional Options
#' @param min_coverage_for_relative If using "relative_cutoff" method, minimum coverage for samples.
#' @param remove_features Logical. Remove features below threshold? Default TRUE.
#' @param generate_plots Logical. Generate QC visualization plots? Default TRUE.
#' @param generate_report Logical. Generate text summary report? Default TRUE.
#' @param verbose Logical. Print progress messages? Default TRUE.
#'
#' @return A list containing:
#'   \item{original_table}{Original loaded table}
#'   \item{filtered_table}{Filtered table}
#'   \item{qc_metrics}{List of QC metrics from compute_filtering_qc()}
#'   \item{presence_stats}{List of presence analysis statistics}
#'   \item{filtering_summary}{Summary of filtering steps applied}
#'
#' @export
#'
#' @examples
#' # Simple pipeline with absolute thresholds
#' # result <- run_filtering_pipeline(
#' #   input_file = "feature-table.tsv",
#' #   output_dir = "results",
#' #   prefix = "analysis1",
#' #   cov_filter_method = "absolute",
#' #   cov_threshold = 1000,
#' #   abun_filter_method = "absolute",
#' #   abun_threshold = 5
#' # )
#'
#' # Pipeline with MAD-based coverage and relative abundance filtering
#' # result <- run_filtering_pipeline(
#' #   input_file = "feature-table.tsv",
#' #   output_dir = "results",
#' #   cov_filter_method = "mad",
#' #   cov_threshold = 3,
#' #   abun_filter_method = "relative_cutoff",
#' #   abun_threshold = 0.01,
#' #   min_coverage_for_relative = 1000
#' # )
run_filtering_pipeline <- function(input_file,
                                    output_dir = ".",
                                    prefix = "filtered",
                                    # Coverage filtering
                                    cov_filter_method = c("none", "absolute", "mad", "iqr", "good", "chao"),
                                    cov_threshold = NULL,
                                    cov_floor = 0,
                                    cov_target_coverage = NULL,
                                    cov_min_reads = 0,
                                    # Abundance filtering
                                    abun_filter_method = c("none", "absolute", "relative", "relative_cutoff"),
                                    abun_threshold = NULL,
                                    abun_min_samples = 1,
                                    # Additional options
                                    min_coverage_for_relative = 1000,
                                    remove_features = TRUE,
                                    generate_plots = TRUE,
                                    generate_report = TRUE,
                                    verbose = TRUE) {
  # Validate inputs
  cov_filter_method <- match.arg(cov_filter_method)
  abun_filter_method <- match.arg(abun_filter_method)

  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file)
  }

  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # === Step 1: Load Table ===
  if (verbose) cat("\n=== Step 1: Loading Table ===\n")
  original_table <- load_feature_table(input_file)
  if (verbose) {
    cat(sprintf("Loaded table with %d features and %d samples\n",
                nrow(original_table), ncol(original_table) - 1))
  }

  # Store original for comparison
  table_current <- original_table

  # Track filtering steps
  filtering_steps <- list()

  # === Step 2: Coverage Filtering ===
  if (cov_filter_method != "none") {
    if (verbose) cat("\n=== Step 2: Coverage Filtering ===\n")

    if (cov_filter_method == "absolute") {
      if (is.null(cov_threshold)) stop("cov_threshold required for absolute method")
      if (verbose) cat(sprintf("Filtering samples with < %d reads\n", cov_threshold))

      table_current <- filter_by_coverage(table_current, min_reads = cov_threshold)
      filtering_steps$coverage <- list(
        method = "absolute",
        threshold = cov_threshold,
        samples_removed = nrow(original_table) - nrow(table_current) + 1  # +1 for feature ID column adjustment
      )

    } else if (cov_filter_method == "mad") {
      if (is.null(cov_threshold)) cov_threshold <- 3
      if (verbose) cat(sprintf("Using MAD method with multiplier = %.2f\n", cov_threshold))

      est <- estimate_mad_cutoff(table_current, multiplier = cov_threshold, floor = cov_floor)
      if (verbose) cat(sprintf("Estimated cutoff: %.0f reads\n", est$cutoff))

      table_current <- filter_by_coverage(table_current, min_reads = est$cutoff)
      filtering_steps$coverage <- list(
        method = "mad",
        multiplier = cov_threshold,
        estimated_cutoff = est$cutoff,
        median = est$median,
        mad = est$mad,
        samples_filtered = est$n_filtered
      )

    } else if (cov_filter_method == "iqr") {
      if (is.null(cov_threshold)) cov_threshold <- 1.5
      if (verbose) cat(sprintf("Using IQR method with multiplier = %.2f\n", cov_threshold))

      est <- estimate_iqr_cutoff(table_current, multiplier = cov_threshold, floor = cov_floor)
      if (verbose) cat(sprintf("Estimated cutoff: %.0f reads\n", est$cutoff))

      table_current <- filter_by_coverage(table_current, min_reads = est$cutoff)
      filtering_steps$coverage <- list(
        method = "iqr",
        multiplier = cov_threshold,
        estimated_cutoff = est$cutoff,
        q1 = est$q1,
        q3 = est$q3,
        iqr = est$iqr,
        samples_filtered = est$n_filtered
      )

    } else if (cov_filter_method == "good") {
      if (is.null(cov_target_coverage)) cov_target_coverage <- 0.95
      if (verbose) cat(sprintf("Using Good's coverage method with target = %.2f (%.0f%%)\n",
                               cov_target_coverage, cov_target_coverage * 100))

      result <- filter_by_coverage_estimator(
        table_current,
        method = "good",
        target_coverage = cov_target_coverage,
        min_reads = cov_min_reads
      )
      table_current <- result$table
      filtering_steps$coverage <- list(
        method = "good",
        target_coverage = cov_target_coverage,
        mean_coverage_before = result$mean_coverage_before,
        mean_coverage_after = result$mean_coverage_after,
        n_samples_before = result$n_samples_before,
        n_samples_after = result$n_samples_after,
        n_samples_filtered = result$n_samples_filtered,
        total_singletons = sum(result$coverage_before == 1 - (1 / colSums(as.matrix(table_current[, -1, drop = FALSE]))), na.rm = TRUE)
      )

    } else if (cov_filter_method == "chao") {
      if (is.null(cov_target_coverage)) cov_target_coverage <- 0.90
      if (verbose) cat(sprintf("Using Chao's coverage method with target = %.2f (%.0f%%)\n",
                               cov_target_coverage, cov_target_coverage * 100))

      result <- filter_by_coverage_estimator(
        table_current,
        method = "chao",
        target_coverage = cov_target_coverage,
        min_reads = cov_min_reads
      )
      table_current <- result$table
      filtering_steps$coverage <- list(
        method = "chao",
        target_coverage = cov_target_coverage,
        mean_coverage_before = result$mean_coverage_before,
        mean_coverage_after = result$mean_coverage_after,
        n_samples_before = result$n_samples_before,
        n_samples_after = result$n_samples_after,
        n_samples_filtered = result$n_samples_filtered
      )
    }

    if (verbose) {
      cat(sprintf("Samples remaining: %d\n", ncol(table_current) - 1))
    }
  } else {
    if (verbose) cat("\n=== Step 2: Coverage Filtering ===\n")
    if (verbose) cat("Skipped (method = none)\n")
    filtering_steps$coverage <- list(method = "none")
  }

  # === Step 3: Abundance Filtering ===
  if (abun_filter_method != "none") {
    if (verbose) cat("\n=== Step 3: Abundance Filtering ===\n")

    if (abun_filter_method == "absolute") {
      if (is.null(abun_threshold)) stop("abun_threshold required for absolute method")
      if (verbose) cat(sprintf("Filtering features with < %d reads\n", abun_threshold))

      table_current <- filter_features_by_abundance(
        table_current,
        threshold = abun_threshold,
        mode = "absolute",
        remove_zeros = remove_features,
        min_samples = abun_min_samples
      )
      filtering_steps$abundance <- list(
        method = "absolute",
        threshold = abun_threshold,
        mode = "absolute",
        features_removed = attr(table_current, "n_filtered_out")
      )

    } else if (abun_filter_method == "relative") {
      if (is.null(abun_threshold)) stop("abun_threshold required for relative method")
      if (verbose) cat(sprintf("Filtering features with < %.4f relative abundance\n", abun_threshold))

      table_current <- filter_features_by_abundance(
        table_current,
        threshold = abun_threshold,
        mode = "relative",
        remove_zeros = remove_features,
        min_samples = abun_min_samples
      )
      filtering_steps$abundance <- list(
        method = "relative",
        threshold = abun_threshold,
        mode = "relative",
        features_removed = attr(table_current, "n_filtered_out")
      )

    } else if (abun_filter_method == "relative_cutoff") {
      if (is.null(abun_threshold)) stop("abun_threshold required for relative_cutoff method")
      if (verbose) cat(sprintf("Using relative cutoff method (%.2f%% of min-coverage sample)\n",
                               abun_threshold * 100))

      result <- filter_by_relative_cutoff(
        table_current,
        min_coverage = min_coverage_for_relative,
        relative_threshold = abun_threshold,
        remove_features = remove_features
      )
      table_current <- result$table
      filtering_steps$abundance <- list(
        method = "relative_cutoff",
        relative_threshold = abun_threshold,
        absolute_threshold = result$absolute_threshold,
        min_sample_coverage = result$min_sample_coverage,
        features_removed = result$n_features_removed
      )
    }

    if (verbose) {
      cat(sprintf("Features remaining: %d\n", nrow(table_current)))
    }
  } else {
    if (verbose) cat("\n=== Step 3: Abundance Filtering ===\n")
    if (verbose) cat("Skipped (method = none)\n")
    filtering_steps$abundance <- list(method = "none")
  }

  # === Step 4: Save Filtered Table ===
  if (verbose) cat("\n=== Step 4: Saving Outputs ===\n")

  # Ensure table is a data.frame with proper column names before writing
  if (is.matrix(table_current)) {
    table_current <- as.data.frame(table_current)
  }
  if (verbose) {
    cat(sprintf("Table type: %s, Columns: %d, First 3 colnames: %s\n",
                class(table_current)[1], ncol(table_current),
                paste(head(colnames(table_current), 3), collapse = ", ")))
  }

  output_file <- file.path(output_dir, paste0(prefix, "_table.tsv"))
  write.table(table_current, file = output_file, sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  if (verbose) cat(sprintf("Filtered table saved to: %s\n", output_file))

  # === Step 5: Generate QC Metrics ===
  if (verbose) cat("\n=== Step 5: Computing QC Metrics ===\n")
  qc_metrics <- compute_filtering_qc(original_table, table_current)

  if (verbose) {
    cat(sprintf("Sparsity drop: %.2f percentage points\n", qc_metrics$sparsity_drop_percent))
    cat(sprintf("Read retention: %.2f%%\n", qc_metrics$read_retention_percent))
    cat(sprintf("Feature retention: %.2f%%\n", qc_metrics$feature_retention_percent))
  }

  # === Step 6: Presence Analysis ===
  if (verbose) cat("\n=== Step 6: Presence Frequency Analysis ===\n")

  # Analyze original vs filtered with comparison plots
  if (verbose) cat("Comparing original vs filtered table...\n")
  presence_stats <- plot_presence_analysis(
    original_table,
    threshold = 1,
    save_dir = if (generate_plots) output_dir else NULL,
    prefix = prefix,
    table_filtered = table_current
  )

  if (verbose) {
    cat(sprintf("Original - Mean sample richness: %.1f features/sample\n", presence_stats$mean_sample_richness))
    cat(sprintf("Filtered - Mean sample richness: %.1f features/sample\n", presence_stats$mean_sample_richness_filtered))
    cat(sprintf("Original - Mean feature prevalence: %.2f samples/feature\n", presence_stats$mean_feature_prevalence))
    cat(sprintf("Filtered - Mean feature prevalence: %.2f samples/feature\n", presence_stats$mean_feature_prevalence_filtered))
  }

  # === Step 7: Generate Plots ===
  if (generate_plots) {
    if (verbose) cat("\n=== Step 7: Generating QC Comparison Plots ===\n")
    plot_paths <- plot_qc_comparison(
      original_table,
      table_current,
      plot_dir = output_dir,
      prefix = prefix
    )
    if (verbose) cat(sprintf("Plots saved to: %s\n", output_dir))
  }

  # === Step 8: Generate Text Report ===
  if (generate_report) {
    if (verbose) cat("\n=== Step 8: Generating Summary Report ===\n")

    report_lines <- c(
      "========================================",
      "FEATURE TABLE FILTERING PIPELINE REPORT",
      "========================================",
      "",
      sprintf("Input file: %s", input_file),
      sprintf("Output directory: %s", output_dir),
      sprintf("Prefix: %s", prefix),
      sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      "",
      "--- ORIGINAL TABLE ---",
      sprintf("Features: %d", nrow(original_table)),
      sprintf("Samples: %d", ncol(original_table) - 1),
      sprintf("Total reads: %s", format(sum(original_table[, -1]), big.mark = ",")),
      "",
      "--- COVERAGE FILTERING ---",
      sprintf("Method: %s", filtering_steps$coverage$method)
    )

    if (filtering_steps$coverage$method == "absolute") {
      report_lines <- c(report_lines,
                        sprintf("Threshold: %d reads", filtering_steps$coverage$threshold))
    } else if (filtering_steps$coverage$method == "mad") {
      report_lines <- c(report_lines,
                        sprintf("Multiplier: %.2f", filtering_steps$coverage$multiplier),
                        sprintf("Median coverage: %.0f", filtering_steps$coverage$median),
                        sprintf("MAD: %.0f", filtering_steps$coverage$mad),
                        sprintf("Cutoff applied: %.0f reads", filtering_steps$coverage$estimated_cutoff))
    } else if (filtering_steps$coverage$method == "iqr") {
      report_lines <- c(report_lines,
                        sprintf("Multiplier: %.2f", filtering_steps$coverage$multiplier),
                        sprintf("Q1: %.0f", filtering_steps$coverage$q1),
                        sprintf("Q3: %.0f", filtering_steps$coverage$q3),
                        sprintf("IQR: %.0f", filtering_steps$coverage$iqr),
                        sprintf("Cutoff applied: %.0f reads", filtering_steps$coverage$estimated_cutoff))
    } else if (filtering_steps$coverage$method == "good") {
      report_lines <- c(report_lines,
                        sprintf("Target coverage: %.2f (%.0f%%)",
                                filtering_steps$coverage$target_coverage,
                                filtering_steps$coverage$target_coverage * 100),
                        sprintf("Mean coverage before filtering: %.4f (%.2f%%)",
                                filtering_steps$coverage$mean_coverage_before,
                                filtering_steps$coverage$mean_coverage_before * 100),
                        sprintf("Mean coverage after filtering: %.4f (%.2f%%)",
                                filtering_steps$coverage$mean_coverage_after,
                                filtering_steps$coverage$mean_coverage_after * 100),
                        sprintf("Samples filtered: %d / %d (%.1f%%)",
                                filtering_steps$coverage$n_samples_filtered,
                                filtering_steps$coverage$n_samples_before,
                                (filtering_steps$coverage$n_samples_filtered / filtering_steps$coverage$n_samples_before) * 100))
    } else if (filtering_steps$coverage$method == "chao") {
      report_lines <- c(report_lines,
                        sprintf("Target coverage: %.2f (%.0f%%)",
                                filtering_steps$coverage$target_coverage,
                                filtering_steps$coverage$target_coverage * 100),
                        sprintf("Mean coverage before filtering: %.4f (%.2f%%)",
                                filtering_steps$coverage$mean_coverage_before,
                                filtering_steps$coverage$mean_coverage_before * 100),
                        sprintf("Mean coverage after filtering: %.4f (%.2f%%)",
                                filtering_steps$coverage$mean_coverage_after,
                                filtering_steps$coverage$mean_coverage_after * 100),
                        sprintf("Samples filtered: %d / %d (%.1f%%)",
                                filtering_steps$coverage$n_samples_filtered,
                                filtering_steps$coverage$n_samples_before,
                                (filtering_steps$coverage$n_samples_filtered / filtering_steps$coverage$n_samples_before) * 100))
    }

    report_lines <- c(report_lines,
                      "",
                      "--- ABUNDANCE FILTERING ---",
                      sprintf("Method: %s", filtering_steps$abundance$method))

    if (filtering_steps$abundance$method == "absolute") {
      report_lines <- c(report_lines,
                        sprintf("Threshold: %d reads", filtering_steps$abundance$threshold))
    } else if (filtering_steps$abundance$method == "relative") {
      report_lines <- c(report_lines,
                        sprintf("Threshold: %.4f (%.2f%%)",
                                filtering_steps$abundance$threshold,
                                filtering_steps$abundance$threshold * 100))
    } else if (filtering_steps$abundance$method == "relative_cutoff") {
      report_lines <- c(report_lines,
                        sprintf("Relative threshold: %.4f (%.2f%%)",
                                filtering_steps$abundance$relative_threshold,
                                filtering_steps$abundance$relative_threshold * 100),
                        sprintf("Absolute threshold applied: %.0f reads",
                                filtering_steps$abundance$absolute_threshold))
    }

    report_lines <- c(report_lines,
                      "",
                      "--- FILTERED TABLE ---",
                      sprintf("Features: %d", nrow(table_current)),
                      sprintf("Samples: %d", ncol(table_current) - 1),
                      sprintf("Total reads: %s", format(sum(table_current[, -1]), big.mark = ",")),
                      "",
                      "--- QC METRICS ---",
                      sprintf("Sparsity (original): %.2f%%", qc_metrics$sparsity_original * 100),
                      sprintf("Sparsity (filtered): %.2f%%", qc_metrics$sparsity_filtered * 100),
                      sprintf("Sparsity drop: %.2f percentage points", qc_metrics$sparsity_drop_percent),
                      sprintf("Read retention: %.2f%%", qc_metrics$read_retention_percent),
                      sprintf("Feature retention: %.2f%%", qc_metrics$feature_retention_percent),
                      sprintf("Sample retention: %.2f%%", qc_metrics$sample_retention_percent),
                      "",
                      "--- PRESENCE ANALYSIS ---",
                      "Original table:",
                      sprintf("  Mean sample richness: %.1f features/sample", presence_stats$mean_sample_richness),
                      sprintf("  Mean feature prevalence: %.2f samples/feature", presence_stats$mean_feature_prevalence),
                      sprintf("  Min/Max richness: %d / %d", presence_stats$min_sample_richness, presence_stats$max_sample_richness),
                      sprintf("  Min/Max prevalence: %d / %d", presence_stats$min_feature_prevalence, presence_stats$max_feature_prevalence),
                      "",
                      "Filtered table:",
                      sprintf("  Mean sample richness: %.1f features/sample", presence_stats$mean_sample_richness_filtered),
                      sprintf("  Mean feature prevalence: %.2f samples/feature", presence_stats$mean_feature_prevalence_filtered),
                      sprintf("  Min/Max richness: %d / %d", presence_stats$min_sample_richness_filtered, presence_stats$max_sample_richness_filtered),
                      sprintf("  Min/Max prevalence: %d / %d", presence_stats$min_feature_prevalence_filtered, presence_stats$max_feature_prevalence_filtered),
                      "",
                      "--- RANK-ABUNDANCE STABILITY ---"
    )

    # Report overlap metrics
    report_lines <- c(report_lines,
                      sprintf("Top N features compared: %d",
                              length(qc_metrics$orig_top_features)),
                      "",
                      "Top features in original table:",
                      paste(sprintf("  %d. %s (%.3f%%)", seq_along(qc_metrics$orig_top_features),
                                    qc_metrics$orig_top_features,
                                    qc_metrics$orig_top_rels * 100)),
                      "",
                      "Top features in filtered table:",
                      paste(sprintf("  %d. %s (%.3f%%)", seq_along(qc_metrics$filt_top_features),
                                    qc_metrics$filt_top_features,
                                    qc_metrics$filt_top_rels * 100)),
                      "",
                      "Overlap between top N lists:",
                      sprintf("  Overlap count: %d / %d",
                              qc_metrics$top_n_overlap_count,
                              max(length(qc_metrics$orig_top_features), length(qc_metrics$filt_top_features))),
                      sprintf("  Overlap percentage: %.1f%%", qc_metrics$top_n_overlap_percent),
                      sprintf("  Jaccard similarity: %.4f", qc_metrics$top_n_jaccard_similarity)
    )

    # Report rank correlation for common features
    if (!is.null(qc_metrics$rank_abundance_correlation) &&
        length(qc_metrics$rank_abundance_correlation) > 0 &&
        !is.na(qc_metrics$rank_abundance_correlation)) {
      report_lines <- c(report_lines,
                        "",
                        "Rank correlation (for features in both top N lists):",
                        sprintf("  Spearman correlation: %.4f", qc_metrics$rank_abundance_correlation),
                        sprintf("  P-value: %.4f", qc_metrics$rank_abundance_pvalue))
    } else {
      report_lines <- c(report_lines,
                        "",
                        "Rank correlation: Not computed (fewer than 2 features in both top N lists)")
    }

    report_lines <- c(report_lines,
                      "",
                      "--- COMPOSITIONAL SIMILARITY (Procrustes) ---"
    )

    if (!is.null(qc_metrics$procrustes_correlation) &&
        length(qc_metrics$procrustes_correlation) > 0 &&
        !is.na(qc_metrics$procrustes_correlation)) {
      report_lines <- c(report_lines,
                        sprintf("Procrustes correlation: %.4f", qc_metrics$procrustes_correlation),
                        sprintf("M^2 statistic: %.6f", qc_metrics$procrustes_m2))
    } else {
      report_lines <- c(report_lines, "Not computed (vegan package not available)")
    }

    report_lines <- c(report_lines,
                      "",
                      "========================================",
                      "END OF REPORT",
                      "========================================"
    )

    report_file <- file.path(output_dir, paste0(prefix, "_report.txt"))
    writeLines(report_lines, report_file)
    if (verbose) cat(sprintf("Report saved to: %s\n", report_file))
  }

  if (verbose) cat("\n=== Pipeline Complete ===\n")

  return(list(
    original_table = original_table,
    filtered_table = table_current,
    qc_metrics = qc_metrics,
    presence_stats = presence_stats,
    filtering_summary = filtering_steps
  ))
}
