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
#' # Singleton Ratio Filtering Options (PCR artifact detection)
#' @param singleton_filter_method Method for singleton ratio filtering:
#'                                "none" - no singleton ratio filtering
#'                                "absolute" - use fixed max_singleton_ratio threshold
#' @param singleton_max_ratio Maximum allowed ratio of singletons+doubletons to total reads.
#'                            Default is 0.1 (10%). Samples exceeding this are removed.
#' @param singleton_count_type Type of low-count features: "singleton", "doubleton", or "both".
#'                             Default is "both".
#'
#' # Cross-Talk / Index Hopping Filtering Options
#' @param crosstalk_filter_method Method for cross-talk filtering:
#'                                "none" - no cross-talk filtering
#'                                "zero" - set suspected leakage reads to zero
#'                                "remove_feature" - remove entire feature if any leakage detected
#'                                "flag" - flag leakage but don't modify data
#' @param crosstalk_threshold Maximum relative abundance threshold for leakage detection.
#'                            Values < this fraction of feature's max are considered leakage.
#'                            Default is 0.001 (0.1% of max).
#' @param crosstalk_min_abs_cutoff Optional minimum absolute count to override relative threshold.
#'                                 Default is NULL (no absolute cutoff).
#' @param crosstalk_return_details Logical. Return detailed leakage matrix? Default FALSE.
#'
#' # Sparsity Elbow Detection Options
#' @param sparsity_elbow_detect Logical. Run sparsity elbow detection for sample filtering recommendation?
#'                              Default TRUE.
#' @param sparsity_elbow_method Method for elbow detection: "kneedle", "max_derivative", or "second_derivative".
#'                              Default is "kneedle".
#' @param apply_sparsity_elbow Logical. Apply sparsity elbow cutoff to filter samples? Default FALSE.
#' @param sparsity_elbow_multiplier MAD multiplier for applying elbow-based coverage cutoff.
#'                                  Default is 1 (use elbow threshold directly).
#'
#' # Depth-Sparsity Outlier Detection Options
#' @param depth_sparsity_detect Logical. Run depth-sparsity outlier analysis? Default TRUE.
#' @param depth_sparsity_metric Metric to analyze: "sparsity" or "richness". Default is "sparsity".
#' @param depth_sparsity_method Outlier detection method: "mad", "iqr", or "both". Default is "mad".
#' @param depth_sparsity_multiplier MAD multiplier for outlier detection. Default is 3.
#' @param depth_sparsity_direction Direction to flag: "high_sparsity", "low_sparsity", or "both".
#'                                  Default is "high_sparsity".
#' @param apply_depth_sparsity_outliers Logical. Remove depth-sparsity outliers? Default FALSE.
#'
#' # Scree Analysis Options (diagnostic)
#' @param run_scree_analysis Logical. Run scree/saturation diagnostic analysis? Default TRUE.
#' @param scree_type Type of scree sweep: "mad_multiplier", "absolute_feature", "relative_feature", or "custom".
#'                   Default is "mad_multiplier".
#' @param scree_n_steps Number of steps for scree analysis. Default is 20.
#' @param scree_custom_thresholds Custom threshold vector for type="custom". Default NULL.
#'
#' # Abundance Filtering Options
#' @param abun_filter_method Method for feature abundance filtering:
#'                           "none" - no abundance filtering
#'                           "absolute" - filter features below absolute count threshold
#'                           "relative" - filter features below relative abundance threshold
#'                           "relative_cutoff" - use relative threshold based on min-coverage sample
#'                           "joint" - joint abundance and prevalence filtering with AND/OR logic
#' @param abun_threshold Threshold value for abundance filtering.
#'                       For "absolute": minimum read count per feature.
#'                       For "relative": proportion (e.g., 0.01 for 1%).
#'                       For "relative_cutoff": relative proportion used to calculate absolute threshold.
#'                       For "joint": minimum abundance threshold for the joint filter.
#' @param abun_min_samples Minimum number of samples where feature must exceed threshold (default: 1).
#'                         Not used for "joint" method (use abun_prevalence_threshold instead).
#'
#' # Joint Filtering Options
#' @param abun_logic Logical operator for joint filtering: "AND" or "OR".
#'                   "AND" requires features to meet BOTH abundance AND prevalence criteria.
#'                   "OR" keeps features meeting EITHER abundance OR prevalence criteria.
#'                   Default is "OR". Only used when abun_filter_method = "joint".
#' @param abun_prevalence_threshold Proportion of samples (0-1) for joint filtering.
#'                                  Features must be present in at least this proportion of samples.
#'                                  Default is 0.3 (30% of samples). Only used when abun_filter_method = "joint".
#'
#' # Additional Options
#' @param min_coverage_for_relative If using "relative_cutoff" method, minimum coverage for samples.
#' @param remove_features Logical. Remove features below threshold? Default TRUE.
#' @param generate_plots Logical. Generate QC visualization plots? Default TRUE.
#' @param generate_report Logical. Generate text summary report? Default TRUE.
#' @param verbose Logical. Print progress messages? Default TRUE.
#'
#' @return A list containing:
#'   \item{original_table}{Original loaded table (same class as input)}
#'   \item{filtered_table}{Filtered table (same class as input)}
#'   \item{qc_metrics}{List of QC metrics from compute_filtering_qc()}
#'   \item{presence_stats}{List of presence analysis statistics}
#'   \item{filtering_summary}{Summary of filtering steps applied}
#'   \item{sparsity_elbow_result}{Result from sparsity elbow detection (if enabled)}
#'   \item{depth_sparsity_result}{Result from depth-sparsity outlier analysis (if enabled)}
#'   \item{scree_result}{Result from scree/saturation analysis (if enabled)}
#'   \item{input_class}{The class of the input object ("data.frame", "phyloseq", or "TreeSummarizedExperiment")}
#'
#' @export
#'
#' @examples
#' # Simple pipeline with absolute thresholds from file
#' # result <- run_filtering_pipeline(
#' #   input = "feature-table.tsv",
#' #   output_dir = "results",
#' #   prefix = "analysis1",
#' #   cov_filter_method = "absolute",
#' #   cov_threshold = 1000,
#' #   abun_filter_method = "absolute",
#' #   abun_threshold = 5
#' # )
#'
#' # Pipeline with phyloseq object (returns phyloseq)
#' # library(phyloseq)
#' # ps <- load_phyloseq("my_data.rds")
#' # result <- run_filtering_pipeline(
#' #   input = ps,
#' #   prefix = "ps_filtered",
#' #   cov_filter_method = "mad",
#' #   cov_threshold = 3,
#' #   abun_filter_method = "relative_cutoff",
#' #   abun_threshold = 0.01
#' # )
#' # filtered_ps <- result$filtered_table  # Still a phyloseq object
#'
#' # Pipeline with TreeSummarizedExperiment (returns TSE)
#' # library(TreeSummarizedExperiment)
#' # tse <- loadTSE("my_data.rds")
#' # result <- run_filtering_pipeline(
#' #   input = tse,
#' #   prefix = "tse_filtered",
#' #   cov_filter_method = "good",
#' #   cov_target_coverage = 0.95,
#' #   abun_filter_method = "joint"
#' # )
#' # filtered_tse <- result$filtered_table  # Still a TreeSummarizedExperiment
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
  abun_filter_method <- match.arg(abun_filter_method)

  # Detect input class and convert to data.frame if necessary
  input_class <- NULL

  if (is.character(input) && length(input) == 1 && file.exists(input)) {
    # Input is a file path
    input_class <- "data.frame"
    if (verbose) cat(sprintf("Loading feature table from file: %s\n", input))
    original_table <- load_feature_table(input)
    input_file_path <- input
  } else if (inherits(input, "phyloseq")) {
    # Input is a phyloseq object
    input_class <- "phyloseq"
    if (verbose) cat("Converting phyloseq object to feature table...\n")
    original_table <- from_phyloseq(input, include_taxa = FALSE)
    input_file_path <- NULL
  } else if (inherits(input, c("TreeSummarizedExperiment", "SingleCellExperiment", "SummarizedExperiment"))) {
    # Input is a TreeSummarizedExperiment or related object
    input_class <- "TreeSummarizedExperiment"
    if (verbose) cat("Converting TreeSummarizedExperiment object to feature table...\n")
    original_table <- from_TSE(input, add_row_data = FALSE)
    input_file_path <- NULL
  } else if (is.data.frame(input) || is.matrix(input)) {
    # Input is already a data.frame or matrix
    input_class <- "data.frame"
    if (verbose) cat("Using provided data.frame/matrix as feature table...\n")
    original_table <- as.data.frame(input)
    if (is.matrix(input)) {
      # Convert matrix to data.frame with feature_id column
      feature_ids <- rownames(original_table)
      if (is.null(feature_ids)) {
        feature_ids <- paste0("Feature_", seq_len(nrow(original_table)))
      }
      original_table <- data.frame(feature_id = feature_ids, original_table, stringsAsFactors = FALSE)
    }
    input_file_path <- NULL
  } else {
    stop("Input must be a file path (character), data.frame, matrix, phyloseq, or TreeSummarizedExperiment object")
  }

  # Create output directory if needed (only for file inputs)
  if (!is.null(input_file_path)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  }

  # === Step 1: Feature Table Ready ===
  if (verbose) cat("\n=== Step 1: Feature Table Prepared ===\n")
  if (verbose) {
    cat(sprintf("Input class: %s\n", input_class))
    cat(sprintf("Table has %d features and %d samples\n",
                nrow(original_table), ncol(original_table) - 1))
  }

  # Store original for comparison
  table_current <- original_table

  # Track filtering steps
  filtering_steps <- list()

  # === Step 2: Singleton Ratio Filtering (PCR artifact detection) ===
  if (singleton_filter_method[1] != "none") {
    if (verbose) cat("\n=== Step 2: Singleton Ratio Filtering ===\n")

    if (singleton_filter_method[1] == "absolute") {
      if (verbose) cat(sprintf("Filtering samples with singleton ratio > %.2f (%s)\n",
                               singleton_max_ratio, singleton_count_type[1]))

      table_current <- filter_by_singleton_ratio(
        table_current,
        max_singleton_ratio = singleton_max_ratio,
        count_type = singleton_count_type[1]
      )

      filtering_steps$singleton_ratio <- list(
        method = singleton_filter_method[1],
        max_ratio = singleton_max_ratio,
        count_type = singleton_count_type[1],
        samples_removed = attr(table_current, "n_filtered_out"),
        samples_retained = attr(table_current, "n_retained")
      )

      if (verbose) {
        cat(sprintf("Samples remaining: %d\n", ncol(table_current) - 1))
      }
    }
  } else {
    if (verbose) cat("\n=== Step 2: Singleton Ratio Filtering ===\n")
    if (verbose) cat("Skipped (method = none)\n")
    filtering_steps$singleton_ratio <- list(method = "none")
  }

  # === Step 3: Coverage Filtering ===
  if (cov_filter_method != "none") {
    if (verbose) cat("\n=== Step 3: Coverage Filtering ===\n")

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
    if (verbose) cat("\n=== Step 3: Coverage Filtering ===\n")
    if (verbose) cat("Skipped (method = none)\n")
    filtering_steps$coverage <- list(method = "none")
  }

  # === Step 4: Cross-Talk / Index Hopping Filtering ===
  if (crosstalk_filter_method[1] != "none") {
    if (verbose) cat("\n=== Step 4: Cross-Talk Filtering ===\n")

    if (verbose) cat(sprintf("Using crosstalk method '%s' with threshold = %.5f\n",
                             crosstalk_filter_method[1], crosstalk_threshold))

    result_ct <- filter_cross_talk(
      table_current,
      max_rel_threshold = crosstalk_threshold,
      min_abs_cutoff = crosstalk_min_abs_cutoff,
      mode = crosstalk_filter_method[1],
      return_details = crosstalk_return_details
    )

    table_current <- result_ct

    filtering_steps$crosstalk <- list(
      method = crosstalk_filter_method[1],
      threshold = crosstalk_threshold,
      min_abs_cutoff = crosstalk_min_abs_cutoff,
      leakage_zeros = attr(result_ct, "n_leakage_zeros"),
      features_affected = attr(result_ct, "n_features_affected")
    )

    if (verbose) {
      cat(sprintf("Features remaining: %d\n", nrow(table_current)))
      if (!is.na(attr(result_ct, "n_leakage_zeros"))) {
        cat(sprintf("Leakage reads zeroed: %d\n", attr(result_ct, "n_leakage_zeros")))
      }
    }
  } else {
    if (verbose) cat("\n=== Step 4: Cross-Talk Filtering ===\n")
    if (verbose) cat("Skipped (method = none)\n")
    filtering_steps$crosstalk <- list(method = "none")
  }

  # === Step 5: Sparsity Elbow Detection ===
  sparsity_elbow_result <- NULL
  if (sparsity_elbow_detect) {
    if (verbose) cat("\n=== Step 5: Sparsity Elbow Detection (Diagnostic) ===\n")

    sparsity_elbow_result <- identify_sparsity_elbow(
      table_current,
      method = sparsity_elbow_method[1]
    )

    if (verbose) {
      cat(sprintf("Elbow threshold: %.0f reads\n", sparsity_elbow_result$elbow_threshold))
      cat(sprintf("Samples above/below: %d / %d\n",
                  sparsity_elbow_result$samples_above_elbow,
                  sparsity_elbow_result$samples_below_elbow))
    }

    filtering_steps$sparsity_elbow <- list(
      detected = TRUE,
      method = sparsity_elbow_method[1],
      elbow_threshold = sparsity_elbow_result$elbow_threshold,
      samples_above = sparsity_elbow_result$samples_above_elbow,
      samples_below = sparsity_elbow_result$samples_below_elbow,
      recommendation = sparsity_elbow_result$recommendation
    )

    if (apply_sparsity_elbow) {
      if (verbose) cat(sprintf("\nApplying sparsity elbow cutoff: %d reads\n",
                               round(sparsity_elbow_result$elbow_threshold * sparsity_elbow_multiplier)))

      table_current <- filter_by_coverage(
        table_current,
        min_reads = sparsity_elbow_result$elbow_threshold * sparsity_elbow_multiplier
      )

      filtering_steps$sparsity_elbow$applied <- TRUE
      filtering_steps$sparsity_elbow$samples_removed_after_apply <-
        sparsity_elbow_result$samples_below_elbow

      if (verbose) {
        cat(sprintf("Samples remaining after elbow filter: %d\n", ncol(table_current) - 1))
      }
    } else {
      filtering_steps$sparsity_elbow$applied <- FALSE
    }
  } else {
    if (verbose) cat("\n=== Step 5: Sparsity Elbow Detection ===\n")
    if (verbose) cat("Skipped (sparsity_elbow_detect = FALSE)\n")
    filtering_steps$sparsity_elbow <- list(detected = FALSE)
  }

  # === Step 6: Depth-Sparsity Outlier Detection ===
  depth_sparsity_result <- NULL
  if (depth_sparsity_detect) {
    if (verbose) cat("\n=== Step 6: Depth-Sparsity Outlier Analysis (Diagnostic) ===\n")

    depth_sparsity_result <- analyze_depth_sparsity(
      table_current,
      metric = depth_sparsity_metric[1],
      outlier_method = depth_sparsity_method[1],
      multiplier = depth_sparsity_multiplier,
      direction = depth_sparsity_direction[1],
      verbose = FALSE  # We'll summarize below
    )

    filtering_steps$depth_sparsity <- list(
      detected = TRUE,
      metric = depth_sparsity_metric[1],
      method = depth_sparsity_method[1],
      n_outliers = depth_sparsity_result$n_outliers,
      outliers = depth_sparsity_result$outliers$sample_name,
      fit_r_squared = depth_sparsity_result$fit_summary$r_squared,
      recommendation = depth_sparsity_result$recommendation
    )

    if (verbose) {
      cat(sprintf("Outliers detected: %d (%.1f%% of samples)\n",
                  depth_sparsity_result$n_outliers,
                  (depth_sparsity_result$n_outliers / (ncol(table_current) - 1)) * 100))
      cat(sprintf("Fit R-squared: %.3f\n", depth_sparsity_result$fit_summary$r_squared))
    }

    if (apply_depth_sparsity_outliers && depth_sparsity_result$n_outliers > 0) {
      if (verbose) cat(sprintf("\nRemoving %d depth-sparsity outliers\n",
                               depth_sparsity_result$n_outliers))

      outlier_names <- depth_sparsity_result$outliers$sample_name
      keep_samples <- setdiff(colnames(table_current)[-1], outlier_names)

      if (length(keep_samples) > 0) {
        # Select feature ID column + kept sample columns by name
        cols_to_keep <- c(colnames(table_current)[1], keep_samples)
        table_current <- table_current[, cols_to_keep, drop = FALSE]

        filtering_steps$depth_sparsity$applied <- TRUE
        filtering_steps$depth_sparsity$samples_removed <- depth_sparsity_result$n_outliers

        if (verbose) {
          cat(sprintf("Samples remaining after outlier removal: %d\n", ncol(table_current) - 1))
        }
      } else {
        if (verbose) cat("Warning: All samples would be removed. Keeping original.\n")
        filtering_steps$depth_sparsity$applied <- FALSE
      }
    } else {
      filtering_steps$depth_sparsity$applied <- FALSE
    }
  } else {
    filtering_steps$depth_sparsity <- list(detected = FALSE)
  }

  # === Step 7: Scree Analysis (Diagnostic) ===
  scree_result <- NULL
  if (run_scree_analysis) {
    if (verbose) cat("\n=== Step 7: Scree/Saturation Analysis ===\n")

    scree_result <- compute_scree(
      table_current,
      type = scree_type[1],
      thresholds = scree_custom_thresholds,
      n_steps = scree_n_steps,
      verbose = verbose
    )

    filtering_steps$scree <- list(
      performed = TRUE,
      type = scree_type[1],
      n_steps = scree_n_steps,
      elbow_threshold = scree_result$summary$elbow_point$threshold,
      elbow_retention = scree_result$summary$elbow_point$retention_at_elbow,
      final_retention = scree_result$summary$saturation$final_retention
    )

    if (verbose) {
      cat(sprintf("Scree analysis complete. Elbow at threshold = %.4f\n",
                  scree_result$summary$elbow_point$threshold))
    }
  } else {
    if (verbose) cat("\n=== Step 7: Scree Analysis ===\n")
    if (verbose) cat("Skipped (run_scree_analysis = FALSE)\n")
    filtering_steps$scree <- list(performed = FALSE)
  }

  # === Step 8: Abundance Filtering ===
  if (abun_filter_method != "none") {
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

    } else if (abun_filter_method == "joint") {
      if (is.null(abun_threshold)) abun_threshold <- 0
      if (verbose) cat(sprintf("Using joint filtering: abundance >= %.4f, prevalence >= %.0f%%, logic = %s\n",
                               abun_threshold, abun_prevalence_threshold * 100, abun_logic))

      result <- filter_features_joint(
        table_current,
        abundance_threshold = abun_threshold,
        prevalence_threshold = abun_prevalence_threshold,
        mode = "relative",
        logic = abun_logic,
        remove_zeros = remove_features
      )
      table_current <- result$table
      filtering_steps$abundance <- list(
        method = "joint",
        abundance_threshold = abun_threshold,
        prevalence_threshold = abun_prevalence_threshold,
        logic = abun_logic,
        n_features_before = result$n_features_before,
        n_features_after = result$n_features_after,
        n_features_removed = result$n_features_removed,
        n_by_abundance_only = result$n_by_abundance_only,
        n_by_prevalence_only = result$n_by_prevalence_only,
        n_by_both = result$n_by_both,
        n_by_neither = result$n_by_neither
      )
    }

    if (verbose) {
      cat(sprintf("Features remaining: %d\n", nrow(table_current)))
    }
  } else {
    if (verbose) cat("\n=== Step 8: Abundance Filtering ===\n")
    if (verbose) cat("Skipped (method = none)\n")
    filtering_steps$abundance <- list(method = "none")
  }

  # === Step 9: Save Filtered Table ===
  if (verbose) cat("\n=== Step 9: Saving Outputs ===\n")

  # Ensure table is a data.frame with proper column names before writing
  if (is.matrix(table_current)) {
    table_current <- as.data.frame(table_current)
  }
  if (verbose) {
    cat(sprintf("Table type: %s, Columns: %d, First 3 colnames: %s\n",
                class(table_current)[1], ncol(table_current),
                paste(head(colnames(table_current), 3), collapse = ", ")))
  }

  # Only save to file if input was a file path
  if (!is.null(input_file_path)) {
    output_file <- file.path(output_dir, paste0(prefix, "_table.tsv"))
    write.table(table_current, file = output_file, sep = "\t",
                row.names = FALSE, col.names = TRUE, quote = FALSE)
    if (verbose) cat(sprintf("Filtered table saved to: %s\n", output_file))
  } else {
    if (verbose) cat("Skipping file output (input was an object)\n")
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

  # === Step 12: Generate Plots ===
  if (generate_plots) {
    if (verbose) cat("\n=== Step 12: Generating QC Comparison Plots ===\n")
    plot_paths <- plot_qc_comparison(
      original_table,
      table_current,
      plot_dir = output_dir,
      prefix = prefix
    )
    if (verbose) cat(sprintf("Plots saved to: %s\n", output_dir))

    # Generate sparsity elbow plot if detection was performed
    if (!is.null(sparsity_elbow_result) && requireNamespace("ggplot2", quietly = TRUE)) {
      if (verbose) cat("Generating Sparsity Elbow plot...\n")
      p_elbow <- plot_sparsity_elbow(sparsity_elbow_result, main = "Sparsity Elbow Detection")
      elbow_plot_path <- file.path(output_dir, paste0(prefix, "_sparsity_elbow.png"))
      ggplot2::ggsave(elbow_plot_path, plot = p_elbow, width = 10, height = 8, dpi = 300)
      if (verbose) cat(sprintf("Sparsity Elbow plot saved to: %s\n", elbow_plot_path))
    }

    # Generate depth-sparsity outlier plot if analysis was performed
    if (!is.null(depth_sparsity_result) && requireNamespace("ggplot2", quietly = TRUE)) {
      if (verbose) cat("Generating Depth-Sparsity Outlier plot...\n")
      p_ds <- plot_depth_sparsity(depth_sparsity_result, main = "Depth-Sparsity Outlier Analysis")
      ds_plot_path <- file.path(output_dir, paste0(prefix, "_depth_sparsity_outliers.png"))
      ggplot2::ggsave(ds_plot_path, plot = p_ds, width = 10, height = 8, dpi = 300)
      if (verbose) cat(sprintf("Depth-Sparsity plot saved to: %s\n", ds_plot_path))
    }

    # Generate scree plot if analysis was performed
    if (!is.null(scree_result) && requireNamespace("ggplot2", quietly = TRUE)) {
      if (verbose) cat("Generating Scree/Saturation plot...\n")
      p_scree <- plot_scree(scree_result, main = "Filtering Threshold Scree Analysis")
      scree_plot_path <- file.path(output_dir, paste0(prefix, "_scree_analysis.png"))
      ggplot2::ggsave(scree_plot_path, plot = p_scree, width = 10, height = 8, dpi = 300)
      if (verbose) cat(sprintf("Scree plot saved to: %s\n", scree_plot_path))
    }
  }

  # === Step 13: Generate Text Report ===
  if (generate_report) {
    if (verbose) cat("\n=== Step 13: Generating Summary Report ===\n")

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
                      "--- SINGLETON RATIO FILTERING ---",
                      sprintf("Method: %s", filtering_steps$singleton_ratio$method))

    if (filtering_steps$singleton_ratio$method == "absolute") {
      report_lines <- c(report_lines,
                        sprintf("Max ratio: %.2f (%.0f%%)",
                                filtering_steps$singleton_ratio$max_ratio,
                                filtering_steps$singleton_ratio$max_ratio * 100),
                        sprintf("Count type: %s", filtering_steps$singleton_ratio$count_type),
                        sprintf("Samples removed: %d", filtering_steps$singleton_ratio$samples_removed))
    }

    report_lines <- c(report_lines,
                      "",
                      "--- CROSS-TALK / INDEX HOPPING FILTERING ---",
                      sprintf("Method: %s", filtering_steps$crosstalk$method))

    if (filtering_steps$crosstalk$method != "none") {
      report_lines <- c(report_lines,
                        sprintf("Threshold: %.5f (%.3f%% of max)", filtering_steps$crosstalk$threshold,
                                filtering_steps$crosstalk$threshold * 100),
                        sprintf("Min absolute cutoff: %s",
                                ifelse(is.null(filtering_steps$crosstalk$min_abs_cutoff), "NULL",
                                       as.character(filtering_steps$crosstalk$min_abs_cutoff))),
                        sprintf("Leakage reads zeroed: %s",
                                ifelse(is.na(filtering_steps$crosstalk$leakage_zeros), "N/A",
                                       as.character(filtering_steps$crosstalk$leakage_zeros))),
                        sprintf("Features affected: %d", filtering_steps$crosstalk$features_affected))
    } else {
      report_lines <- c(report_lines, "No cross-talk filtering applied")
    }

    report_lines <- c(report_lines,
                      "",
                      "--- SPARSITY ELBOW DETECTION ---",
                      sprintf("Detection performed: %s",
                              ifelse(filtering_steps$sparsity_elbow$detected, "Yes", "No")))

    if (filtering_steps$sparsity_elbow$detected) {
      report_lines <- c(report_lines,
                        sprintf("Method: %s", filtering_steps$sparsity_elbow$method),
                        sprintf("Elbow threshold: %.0f reads", filtering_steps$sparsity_elbow$elbow_threshold),
                        sprintf("Samples above elbow: %d", filtering_steps$sparsity_elbow$samples_above),
                        sprintf("Samples below elbow: %d", filtering_steps$sparsity_elbow$samples_below),
                        sprintf("Applied to filter: %s",
                                ifelse(filtering_steps$sparsity_elbow$applied, "Yes", "No")))
    }

    report_lines <- c(report_lines,
                      "",
                      "--- DEPTH-SPARSITY OUTLIER ANALYSIS ---",
                      sprintf("Analysis performed: %s",
                              ifelse(filtering_steps$depth_sparsity$detected, "Yes", "No")))

    if (filtering_steps$depth_sparsity$detected) {
      report_lines <- c(report_lines,
                        sprintf("Metric: %s", filtering_steps$depth_sparsity$metric),
                        sprintf("Method: %s", filtering_steps$depth_sparsity$method),
                        sprintf("Outliers detected: %d", filtering_steps$depth_sparsity$n_outliers),
                        sprintf("Fit R-squared: %.3f", filtering_steps$depth_sparsity$fit_r_squared),
                        sprintf("Applied to filter: %s",
                                ifelse(!is.null(filtering_steps$depth_sparsity$applied) &&
                                       filtering_steps$depth_sparsity$applied, "Yes", "No")))

      if (filtering_steps$depth_sparsity$n_outliers > 0) {
        outlier_list <- paste(head(filtering_steps$depth_sparsity$outliers, 10), collapse = ", ")
        if (length(filtering_steps$depth_sparsity$outliers) > 10) {
          outlier_list <- paste0(outlier_list, " ... (and more)")
        }
        report_lines <- c(report_lines,
                          sprintf("Outlier samples: %s", outlier_list))
      }
    }

    report_lines <- c(report_lines,
                      "",
                      "--- SCREE / SATURATION ANALYSIS ---",
                      sprintf("Analysis performed: %s",
                              ifelse(filtering_steps$scree$performed, "Yes", "No")))

    if (filtering_steps$scree$performed) {
      report_lines <- c(report_lines,
                        sprintf("Type: %s", filtering_steps$scree$type),
                        sprintf("Steps evaluated: %d", filtering_steps$scree$n_steps),
                        sprintf("Elbow threshold: %.4f", filtering_steps$scree$elbow_threshold),
                        sprintf("Retention at elbow: %.1f%%", filtering_steps$scree$elbow_retention),
                        sprintf("Final retention: %.1f%%", filtering_steps$scree$final_retention))
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
    } else if (filtering_steps$abundance$method == "joint") {
      report_lines <- c(report_lines,
                        sprintf("Abundance threshold: %.4f (%.2f%%)",
                                filtering_steps$abundance$abundance_threshold,
                                filtering_steps$abundance$abundance_threshold * 100),
                        sprintf("Prevalence threshold: %.2f (%.0f%% of samples)",
                                filtering_steps$abundance$prevalence_threshold,
                                filtering_steps$abundance$prevalence_threshold * 100),
                        sprintf("Logic: %s", filtering_steps$abundance$logic),
                        "",
                        "Feature retention breakdown:",
                        sprintf("  Met abundance only: %d features", filtering_steps$abundance$n_by_abundance_only),
                        sprintf("  Met prevalence only: %d features", filtering_steps$abundance$n_by_prevalence_only),
                        sprintf("  Met both criteria: %d features", filtering_steps$abundance$n_by_both),
                        sprintf("  Met neither (removed): %d features", filtering_steps$abundance$n_by_neither),
                        "",
                        sprintf("Features before: %d", filtering_steps$abundance$n_features_before),
                        sprintf("Features after: %d", filtering_steps$abundance$n_features_after),
                        sprintf("Features removed: %d (%.1f%%)",
                                filtering_steps$abundance$n_features_removed,
                                (filtering_steps$abundance$n_features_removed / filtering_steps$abundance$n_features_before) * 100))
    }

    report_lines <- c(report_lines,
                      "",
                      "--- FILTERED TABLE ---",
                      sprintf("Features: %d", nrow(table_current)),
                      sprintf("Samples: %d", ncol(table_current) - 1),
                      sprintf("Total reads: %s", format(as.integer(sum(table_current[, -1])), big.mark = ",")),
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

    # Add diagnostic plots section if generated
    report_lines <- c(report_lines,
                      "",
                      "--- DIAGNOSTIC PLOTS ---"
    )

    if (generate_plots) {
      if (!is.null(sparsity_elbow_result)) {
        report_lines <- c(report_lines,
                          sprintf("Sparsity Elbow plot: %s", paste0(prefix, "_sparsity_elbow.png")))
      }
      if (!is.null(depth_sparsity_result)) {
        report_lines <- c(report_lines,
                          sprintf("Depth-Sparsity Outlier plot: %s", paste0(prefix, "_depth_sparsity_outliers.png")))
      }
      if (!is.null(scree_result)) {
        report_lines <- c(report_lines,
                          sprintf("Scree/Saturation plot: %s", paste0(prefix, "_scree_analysis.png")))
      }
    } else {
      report_lines <- c(report_lines, "Plots were not generated (generate_plots = FALSE)")
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

  # === Convert Output Back to Original Input Class ===
  if (verbose) cat(sprintf("\nConverting output back to %s format...\n", input_class))

  # Convert original and filtered tables back to original input class
  if (input_class == "phyloseq") {
    # Get any additional components from the original phyloseq object
    original_ps <- input
    filtered_ps <- to_phyloseq(table_current,
                                tax_table = phyloseq::tax_table(original_ps),
                                phy_tree = phyloseq::phy_tree(original_ps),
                                sample_data = phyloseq::sample_data(original_ps))
    original_table <- original_ps
    table_current <- filtered_ps
  } else if (input_class == "TreeSummarizedExperiment") {
    # Get any additional components from the original TSE object
    original_tse <- input
    filtered_tse <- to_TSE(table_current,
                            rowData = SummarizedExperiment::rowData(original_tse),
                            colData = SummarizedExperiment::colData(original_tse),
                            reducedDims = TreeSummarizedExperiment::reducedDims(original_tse),
                            rowTree = TreeSummarizedExperiment::rowTree(original_tse),
                            rowLinks = TreeSummarizedExperiment::rowLinks(original_tse))
    original_table <- original_tse
    table_current <- filtered_tse
  }
  # For data.frame input, no conversion needed

  if (verbose) cat("\n=== Pipeline Complete ===\n")

  return(list(
    original_table = original_table,
    filtered_table = table_current,
    qc_metrics = qc_metrics,
    presence_stats = presence_stats,
    filtering_summary = filtering_steps,
    sparsity_elbow_result = sparsity_elbow_result,
    depth_sparsity_result = depth_sparsity_result,
    scree_result = scree_result,
    input_class = input_class
  ))
}
