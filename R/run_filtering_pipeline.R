#' Run a complete feature table filtering pipeline with QC report
#'
#' A comprehensive pipeline that loads, filters, and analyzes a feature table.
#' Allows selection of different coverage, abundance, and quality-based filtering methods.
#' Produces QC metrics and visualizations comparing before/after filtering.
#' Supports data.frame, phyloseq, and TreeSummarizedExperiment inputs and returns
#' the same class as the input.
#'
#' @param input Path to input feature table file (TSV or CSV), or a feature table object
#'              (data.frame, phyloseq, or TreeSummarizedExperiment).
#' @param output_dir Directory to save outputs (filtered table, plots, report). Default is current directory.
#'                   Ignored when input is an object rather than a file path.
#' @param prefix Prefix for output filenames. Default is "filtered".
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
#' @param singleton_filter_method Method for singleton ratio filtering:
#'                                "none" - no singleton ratio filtering
#'                                "absolute" - use fixed max_singleton_ratio threshold
#' @param singleton_max_ratio Maximum allowed ratio of singletons+doubletons to total reads.
#'                            Default is 0.1 (10\%). Samples exceeding this are removed.
#' @param singleton_count_type Type of low-count features: "singleton", "doubleton", or "both".
#'                             Default is "both".
#' @param crosstalk_filter_method Method for cross-talk filtering:
#'                                "none" - no cross-talk filtering
#'                                "zero" - set suspected leakage reads to zero
#'                                "remove_feature" - remove entire feature if any leakage detected
#'                                "flag" - flag leakage but don't modify data
#' @param crosstalk_threshold Maximum relative abundance threshold for leakage detection.
#'                            Values < this fraction of feature's max are considered leakage.
#'                            Default is 0.001 (0.1\% of max).
#' @param crosstalk_min_abs_cutoff Optional minimum absolute count to override relative threshold.
#'                                 Default is NULL (no absolute cutoff).
#' @param crosstalk_return_details Logical. Return detailed leakage matrix? Default FALSE.
#' @param sparsity_elbow_detect Logical. Run sparsity elbow detection for sample filtering recommendation?
#'                              Default TRUE.
#' @param sparsity_elbow_method Method for elbow detection: "kneedle", "max_derivative", or "second_derivative".
#'                              Default is "kneedle".
#' @param apply_sparsity_elbow Logical. Apply sparsity elbow cutoff to filter samples? Default FALSE.
#' @param sparsity_elbow_multiplier MAD multiplier for applying elbow-based coverage cutoff.
#'                                  Default is 1 (use elbow threshold directly).
#' @param depth_sparsity_detect Logical. Run depth-sparsity outlier analysis? Default TRUE.
#' @param depth_sparsity_metric Metric to analyze: "sparsity" or "richness". Default is "sparsity".
#' @param depth_sparsity_method Outlier detection method: "mad", "iqr", or "both". Default is "mad".
#' @param depth_sparsity_multiplier MAD multiplier for outlier detection. Default is 3.
#' @param depth_sparsity_direction Direction to flag: "high_sparsity", "low_sparsity", or "both".
#'                                  Default is "high_sparsity".
#' @param apply_depth_sparsity_outliers Logical. Remove depth-sparsity outliers? Default FALSE.
#' @param run_scree_analysis Logical. Run scree/saturation diagnostic analysis? Default TRUE.
#' @param scree_type Type of scree sweep: "mad_multiplier", "absolute_feature", "relative_feature", or "custom".
#'                   Default is "mad_multiplier".
#' @param scree_n_steps Number of steps for scree analysis. Default is 20.
#' @param scree_custom_thresholds Custom threshold vector for type="custom". Default NULL.
#' @param abun_filter_method Method for feature abundance filtering:
#'                           "none" - no abundance filtering
#'                           "absolute" - filter features below absolute count threshold
#'                           "relative" - filter features below relative abundance threshold
#'                           "relative_cutoff" - use relative threshold based on min-coverage sample
#'                           "joint" - joint abundance and prevalence filtering with AND/OR logic
#' @param abun_threshold Threshold value for abundance filtering.
#'                       For "absolute": minimum read count per feature.
#'                       For "relative": proportion (e.g., 0.01 for 1\%).
#'                       For "relative_cutoff": relative proportion used to calculate absolute threshold.
#'                       For "joint": minimum abundance threshold for the joint filter.
#' @param abun_min_samples Minimum number of samples where feature must exceed threshold (default: 1).
#'                         Not used for "joint" method (use abun_prevalence_threshold instead).
#' @param abun_logic Logical operator for joint filtering: "AND" or "OR".
#'                   "AND" requires features to meet BOTH abundance AND prevalence criteria.
#'                   "OR" keeps features meeting EITHER abundance OR prevalence criteria.
#'                   Default is "OR". Only used when abun_filter_method = "joint".
#' @param abun_prevalence_threshold Proportion of samples (0-1) for joint filtering.
#'                                  Features must be present in at least this proportion of samples.
#'                                  Default is 0.3 (30\% of samples). Only used when abun_filter_method = "joint".
#' @param min_coverage_for_relative If using "relative_cutoff" method, minimum coverage for samples.
#' @param remove_features Logical. Remove features below threshold? Default TRUE.
#' @param generate_plots Logical. Generate QC visualization plots? Default TRUE.
#' @param generate_report Logical. Generate text summary report? Default TRUE.
#' @param verbose Logical. Print progress messages? Default TRUE.
#'
#' @return A list containing:
#' \describe{
#'   \item{original_table}{Original loaded table (same class as input)}
#'   \item{filtered_table}{Filtered table (same class as input)}
#'   \item{qc_metrics}{List of QC metrics from compute_filtering_qc()}
#'   \item{presence_stats}{List of presence analysis statistics}
#'   \item{filtering_summary}{Summary of filtering steps applied}
#'   \item{sparsity_elbow_result}{Result from sparsity elbow detection (if enabled)}
#'   \item{depth_sparsity_result}{Result from depth-sparsity outlier analysis (if enabled)}
#'   \item{scree_result}{Result from scree/saturation analysis (if enabled)}
#'   \item{input_class}{The class of the input object ("data.frame", "phyloseq", or "TreeSummarizedExperiment")}
#' }
#'
#' @export
#'
#' @examples
#' data(example_feature_table)
#' result <- run_filtering_pipeline(example_feature_table,
#'                                  cov_filter_method = "absolute",
#'                                  cov_threshold = 1000,
#'                                  abun_filter_method = "absolute",
#'                                  abun_threshold = 10,
#'                                  generate_report = FALSE,
#'                                  verbose = FALSE)
#' nrow(result$filtered_table)
run_filtering_pipeline <- function(input,
                                    output_dir = ".",
                                    prefix = "filtered",
                                    # Coverage filtering
                                    cov_filter_method = c("none", "absolute", "mad", "iqr", "good", "chao"),
                                    cov_threshold = NULL,
                                    cov_floor = 0,
                                    cov_target_coverage = NULL,
                                    cov_min_reads = 0,
                                    # Singleton ratio filtering
                                    singleton_filter_method = c("none", "absolute"),
                                    singleton_max_ratio = 0.1,
                                    singleton_count_type = c("both", "singleton", "doubleton"),
                                    # Cross-talk filtering
                                    crosstalk_filter_method = c("none", "zero", "remove_feature", "flag"),
                                    crosstalk_threshold = 0.001,
                                    crosstalk_min_abs_cutoff = NULL,
                                    crosstalk_return_details = FALSE,
                                    # Sparsity elbow detection
                                    sparsity_elbow_detect = TRUE,
                                    sparsity_elbow_method = c("kneedle", "max_derivative", "second_derivative"),
                                    apply_sparsity_elbow = FALSE,
                                    sparsity_elbow_multiplier = 1,
                                    # Depth-sparsity outlier detection
                                    depth_sparsity_detect = TRUE,
                                    depth_sparsity_metric = c("sparsity", "richness"),
                                    depth_sparsity_method = c("mad", "iqr", "both"),
                                    depth_sparsity_multiplier = 3,
                                    depth_sparsity_direction = c("high_sparsity", "low_sparsity", "both"),
                                    apply_depth_sparsity_outliers = FALSE,
                                    # Scree analysis
                                    run_scree_analysis = TRUE,
                                    scree_type = c("mad_multiplier", "absolute_feature", "relative_feature", "custom"),
                                    scree_n_steps = 20,
                                    scree_custom_thresholds = NULL,
                                    # Abundance filtering
                                    abun_filter_method = c("none", "absolute", "relative", "relative_cutoff", "joint"),
                                    abun_threshold = NULL,
                                    abun_min_samples = 1,
                                    # Joint filtering options
                                    abun_logic = c("OR", "AND"),
                                    abun_prevalence_threshold = 0.3,
                                    # Additional options
                                    min_coverage_for_relative = 1000,
                                    remove_features = TRUE,
                                    generate_plots = TRUE,
                                    generate_report = TRUE,
                                    verbose = TRUE) {
  # Validate inputs
  cov_filter_method <- match.arg(cov_filter_method)
  singleton_filter_method <- match.arg(singleton_filter_method)
  singleton_count_type <- match.arg(singleton_count_type)
  crosstalk_filter_method <- match.arg(crosstalk_filter_method)
  sparsity_elbow_method <- match.arg(sparsity_elbow_method)
  depth_sparsity_metric <- match.arg(depth_sparsity_metric)
  depth_sparsity_method <- match.arg(depth_sparsity_method)
  depth_sparsity_direction <- match.arg(depth_sparsity_direction)
  scree_type <- match.arg(scree_type)
  abun_filter_method <- match.arg(abun_filter_method)
  abun_logic <- match.arg(abun_logic)

  # === Step 1: Load Input Data ===
  if (verbose) cat("\n=== Step 1: Loading Input Data ===\n")
  input_data <- .load_pipeline_input(input, verbose = verbose)
  original_table <- input_data$table
  input_class <- input_data$input_class
  input_file_path <- input_data$input_file_path
  original_object <- input_data$original_object

  # Store original stats
  original_stats <- .get_table_stats(original_table)
  table_current <- original_table

  # Create output directory if needed
  output_dir <- .create_output_directory(output_dir, input_file_path, verbose)

  # Track filtering steps
  filtering_steps <- list()

  # === Step 2: Singleton Ratio Filtering ===
  if (singleton_filter_method != "none") {
    if (verbose) cat("\n=== Step 2: Singleton Ratio Filtering ===\n")
    stats_before <- .get_table_stats(table_current)
    table_current <- .apply_singleton_filter(
      table_current, singleton_max_ratio, singleton_count_type, verbose
    )
    stats_after <- .get_table_stats(table_current)
    filtering_steps$singleton_ratio <- .format_step_summary(
      "Singleton Ratio", singleton_filter_method,
      list(max_ratio = singleton_max_ratio, count_type = singleton_count_type),
      stats_before, stats_after
    )
  }

  # === Step 3: Coverage Filtering ===
  if (cov_filter_method != "none") {
    if (verbose) cat("\n=== Step 3: Coverage Filtering ===\n")
    stats_before <- .get_table_stats(table_current)
    table_current <- .apply_coverage_filter(
      table_current, cov_filter_method, cov_threshold, cov_floor,
      cov_target_coverage, cov_min_reads, verbose
    )
    stats_after <- .get_table_stats(table_current)
    filtering_steps$coverage <- .format_step_summary(
      "Coverage", cov_filter_method,
      list(threshold = cov_threshold, floor = cov_floor),
      stats_before, stats_after
    )
  }

  # === Step 4: Cross-Talk Filtering ===
  if (crosstalk_filter_method != "none") {
    if (verbose) cat("\n=== Step 4: Cross-Talk Filtering ===\n")
    stats_before <- .get_table_stats(table_current)
    crosstalk_result <- .apply_crosstalk_filter(
      table_current, crosstalk_filter_method, crosstalk_threshold,
      crosstalk_min_abs_cutoff, crosstalk_return_details, verbose
    )
    stats_after <- .get_table_stats(crosstalk_result$table)
    table_current <- crosstalk_result$table
    filtering_steps$crosstalk <- .format_step_summary(
      "Cross-Talk", crosstalk_filter_method,
      list(threshold = crosstalk_threshold, min_abs_cutoff = crosstalk_min_abs_cutoff),
      stats_before, stats_after
    )
    # Store detailed leakage info if returned
    if (!is.null(crosstalk_result$detailed_leakage)) {
      filtering_steps$crosstalk$detailed_leakage <- crosstalk_result$detailed_leakage
    }
  }

  # === Step 5: Sparsity Elbow Detection/Filtering ===
  sparsity_elbow_result <- NULL
  if (sparsity_elbow_detect || apply_sparsity_elbow) {
    if (verbose) cat("\n=== Step 5: Sparsity Elbow Analysis ===\n")
    elbow_result <- .run_sparsity_elbow_analysis(
      table_current, sparsity_elbow_method, verbose
    )
    sparsity_elbow_result <- elbow_result

    if (apply_sparsity_elbow) {
      stats_before <- .get_table_stats(table_current)
      # Apply elbow-based coverage filter
      elbow_threshold <- elbow_result$kneedle_threshold
      if (!is.null(elbow_threshold) && !is.na(elbow_threshold)) {
        adjusted_cutoff <- max(elbow_threshold - sparsity_elbow_multiplier * mad(elbow_result$richness_values), 0)
        table_current <- filter_by_coverage(table_current, min_reads = adjusted_cutoff)
      }
      stats_after <- .get_table_stats(table_current)
      filtering_steps$sparsity_elbow <- .format_step_summary(
        "Sparsity Elbow", sparsity_elbow_method,
        list(multiplier = sparsity_elbow_multiplier),
        stats_before, stats_after
      )
    }
  }

  # === Step 6: Depth-Sparsity Outlier Detection/Filtering ===
  depth_sparsity_result <- NULL
  if (depth_sparsity_detect || apply_depth_sparsity_outliers) {
    if (verbose) cat("\n=== Step 6: Depth-Sparsity Analysis ===\n")
    ds_result <- .run_depth_sparsity_analysis(
      table_current, depth_sparsity_metric, depth_sparsity_method,
      depth_sparsity_multiplier, depth_sparsity_direction, verbose
    )
    depth_sparsity_result <- ds_result

    if (apply_depth_sparsity_outliers) {
      stats_before <- .get_table_stats(table_current)
      table_current <- filter_depth_sparsity_outliers(
        table_current,
        metric = depth_sparsity_metric,
        outlier_method = depth_sparsity_method,
        multiplier = depth_sparsity_multiplier,
        direction = depth_sparsity_direction
      )
      stats_after <- .get_table_stats(table_current)
      filtering_steps$depth_sparsity <- .format_step_summary(
        "Depth-Sparsity Outliers", depth_sparsity_method,
        list(metric = depth_sparsity_metric, multiplier = depth_sparsity_multiplier),
        stats_before, stats_after
      )
    }
  }

  # === Step 7: Scree Analysis ===
  scree_result <- NULL
  if (run_scree_analysis) {
    if (verbose) cat("\n=== Step 7: Scree Analysis ===\n")
    scree_result <- .run_scree_analysis(
      original_table, scree_type, scree_n_steps,
      scree_custom_thresholds, verbose
    )
  }

  # === Step 8: Abundance Filtering ===
  if (abun_filter_method != "none") {
    if (verbose) cat("\n=== Step 8: Abundance Filtering ===\n")
    stats_before <- .get_table_stats(table_current)
    table_current <- .apply_abundance_filter(
      table_current, abun_filter_method, abun_threshold, abun_min_samples,
      abun_logic, abun_prevalence_threshold, min_coverage_for_relative,
      remove_features, verbose
    )
    stats_after <- .get_table_stats(table_current)
    filtering_steps$abundance <- .format_step_summary(
      "Abundance", abun_filter_method,
      list(threshold = abun_threshold, min_samples = abun_min_samples),
      stats_before, stats_after
    )
  }

  # === Step 8b: Remove any remaining all-zero features ===
  stats_before <- .get_table_stats(table_current)
  sample_cols <- as.matrix(table_current[, -1, drop = FALSE])
  keep_features <- rowSums(sample_cols) > 0
  n_removed <- sum(!keep_features)
  if (n_removed > 0) {
    if (verbose) cat(sprintf("Removing %d features with all zeros...\n", n_removed))
    table_current <- table_current[keep_features, , drop = FALSE]
  }
  stats_after <- .get_table_stats(table_current)

  # === Step 9: Save Filtered Table ===
  if (verbose) cat("\n=== Step 9: Saving Outputs ===\n")
  if (!is.null(input_file_path)) {
    output_file <- file.path(output_dir, paste0(prefix, "_table.tsv"))
    write.table(table_current, file = output_file, sep = "\t",
                row.names = FALSE, col.names = TRUE, quote = FALSE)
    if (verbose) cat(sprintf("Filtered table saved to: %s\n", output_file))
  }

  # === Step 10: Compute QC Metrics ===
  if (verbose) cat("\n=== Step 10: Computing QC Metrics ===\n")
  qc_metrics <- compute_filtering_qc(original_table, table_current)
  if (verbose) {
    cat(sprintf("Sparsity drop: %.2f percentage points\n", qc_metrics$sparsity_drop_percent))
    cat(sprintf("Read retention: %.2f%%\n", qc_metrics$read_retention_percent))
    cat(sprintf("Feature retention: %.2f%%\n", qc_metrics$feature_retention_percent))
  }

  # === Step 11: Presence Analysis ===
  if (verbose) cat("\n=== Step 11: Presence Frequency Analysis ===\n")
  presence_stats <- plot_presence_analysis(
    original_table, threshold = 1,
    save_dir = if (generate_plots) output_dir else NULL,
    prefix = prefix, table_filtered = table_current
  )

  # === Step 12: Generate Plots ===
  if (generate_plots) {
    if (verbose) cat("\n=== Step 12: Generating Plots ===\n")
    plot_paths <- .generate_qc_plots(original_table, table_current, output_dir, prefix, verbose)

    if (!is.null(sparsity_elbow_result) && requireNamespace("ggplot2", quietly = TRUE)) {
      .generate_sparsity_elbow_plot(sparsity_elbow_result, output_dir, prefix,
                                     "Sparsity Elbow Detection", verbose)
    }
    if (!is.null(depth_sparsity_result) && requireNamespace("ggplot2", quietly = TRUE)) {
      .generate_depth_sparsity_plot(depth_sparsity_result, output_dir, prefix,
                                     "Depth-Sparsity Outlier Analysis", verbose)
    }
    if (!is.null(scree_result) && requireNamespace("ggplot2", quietly = TRUE)) {
      .generate_scree_plot(scree_result, output_dir, prefix,
                           "Filtering Threshold Scree Analysis", verbose)
    }
  }

  # === Step 13: Generate Report ===
  if (generate_report) {
    if (verbose) cat("\n=== Step 13: Generating Report ===\n")
    # Determine input description
    if (!is.null(input_file_path)) {
      input_desc <- input_file_path
    } else {
      input_desc <- sprintf("%s object (%d features x %d samples)",
                            input_class, original_stats$features, original_stats$samples)
    }
    # Collect pipeline parameters for report
    pipeline_params <- list(
      cov_filter_method = cov_filter_method,
      cov_threshold = cov_threshold,
      cov_floor = cov_floor,
      cov_target_coverage = cov_target_coverage,
      cov_min_reads = cov_min_reads,
      singleton_filter_method = singleton_filter_method,
      singleton_max_ratio = singleton_max_ratio,
      crosstalk_filter_method = crosstalk_filter_method,
      crosstalk_threshold = crosstalk_threshold,
      apply_sparsity_elbow = apply_sparsity_elbow,
      apply_depth_sparsity_outliers = apply_depth_sparsity_outliers,
      abun_filter_method = abun_filter_method,
      abun_threshold = abun_threshold,
      abun_min_samples = abun_min_samples,
      abun_logic = abun_logic,
      abun_prevalence_threshold = abun_prevalence_threshold
    )
    # Get filtered table stats
    filtered_stats <- .get_table_stats(table_current)

    report_paths <- .generate_filtering_report(
      original_stats, qc_metrics, presence_stats, filtering_steps,
      sparsity_elbow_result, depth_sparsity_result, scree_result,
      output_dir, prefix, verbose,
      input_description = input_desc,
      pipeline_params = pipeline_params,
      filtered_stats = filtered_stats
    )
  }

  # === Convert Output Back to Original Class ===
  if (verbose) cat(sprintf("\nConverting output back to %s format...\n", input_class))
  converted <- .convert_output_back_to_original_class(
    table_current, input_class, original_object, original_table
  )

  if (verbose) cat("\n=== Pipeline Complete ===\n")

  return(list(
    original_table = converted$original_table,
    filtered_table = converted$filtered_table,
    qc_metrics = qc_metrics,
    presence_stats = presence_stats,
    filtering_summary = filtering_steps,
    sparsity_elbow_result = sparsity_elbow_result,
    depth_sparsity_result = depth_sparsity_result,
    scree_result = scree_result,
    input_class = input_class
  ))
}
