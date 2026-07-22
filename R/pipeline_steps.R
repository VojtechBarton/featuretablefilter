# =============================================================================
# Pipeline Filtering Step Wrapper Functions (Internal)
# =============================================================================
# These functions wrap the existing filtering functions with standardized
# interfaces for use in the pipeline. All functions are internal (. prefix).
# =============================================================================

#' Apply singleton ratio filtering#'#' Wrapper around filter_by_singleton_ratio() with standardized return format.#'#' @param table Feature table (data.frame)#' @param max_ratio Maximum allowed singleton ratio#' @param count_type Type of low-count features: "singleton", "doubleton", or "both"#' @param verbose Logical. Print progress messages?#'#' @return Filtered feature table
#' @noRd
.apply_singleton_filter <- function(table, max_ratio, count_type, verbose = TRUE) {
  if (verbose) cat(sprintf("Applying singleton ratio filter (max ratio: %.2f, type: %s)...\n",
                           max_ratio, count_type))
  filter_by_singleton_ratio(table, max_singleton_ratio = max_ratio, count_type = count_type)
}

#' Apply coverage filtering#'#' Wrapper that handles different coverage filtering methods and returns#' a standardized result.#'#' @param table Feature table (data.frame)#' @param method Method: "absolute", "mad", "iqr", "good", "chao"#' @param threshold Fixed threshold or multiplier#' @param floor Minimum possible cutoff#' @param target_coverage Target ecological coverage (for good/chao methods)#' @param min_reads Optional minimum absolute read count floor#' @param verbose Logical. Print progress messages?#'#' @return Filtered feature table
#' @noRd
.apply_coverage_filter <- function(table, method, threshold = NULL, floor = 0,
                                    target_coverage = NULL, min_reads = 0, verbose = TRUE) {
  if (method == "none") {
    if (verbose) cat("Skipping coverage filtering (method = 'none')\n")
    return(table)
  }

  # Determine cutoff value based on method
  cutoff <- switch(method,
    "absolute" = threshold,
    "mad" = {
      est <- estimate_mad_cutoff(table, multiplier = threshold, floor = floor)
      max(est$cutoff, min_reads)
    },
    "iqr" = {
      est <- estimate_iqr_cutoff(table, multiplier = threshold, floor = floor)
      max(est$cutoff, min_reads)
    },
    "good" = {
      est <- filter_by_coverage_estimator(table, method = "good",
                                           target_coverage = target_coverage,
                                           min_reads = min_reads,
                                           verbose = verbose)
      return(est$table)
    },
    "chao" = {
      est <- filter_by_coverage_estimator(table, method = "chao",
                                           target_coverage = target_coverage,
                                           min_reads = min_reads,
                                           verbose = verbose)
      return(est$table)
    },
    stop("Unknown coverage method: ", method)
  )

  if (verbose) cat(sprintf("Filtering samples with coverage < %.0f reads...\n", cutoff))
  filter_by_coverage(table, min_reads = cutoff)
}

#' Apply cross-talk (index hopping) filtering#'#' Wrapper around filter_cross_talk() with standardized interface.#'#' @param table Feature table (data.frame)#' @param method Method: "zero", "remove_feature", "flag", or "none"#' @param threshold Maximum relative abundance threshold#' @param min_abs_cutoff Minimum absolute count override#' @param return_details Return detailed leakage matrix?#' @param verbose Logical. Print progress messages?#'#' @return List with filtered table and optional detailed leakage info
#' @noRd
.apply_crosstalk_filter <- function(table, method, threshold, min_abs_cutoff,
                                     return_details, verbose = TRUE) {
  if (method == "none") {
    if (verbose) cat("Skipping cross-talk filtering (method = 'none')\n")
    return(list(table = table, detailed_leakage = NULL))
  }

  if (verbose) cat(sprintf("Applying cross-talk filter (threshold: %.4f, method: %s)...\n",
                           threshold, method))
  result <- filter_cross_talk(table, max_rel_threshold = threshold,
                              min_abs_cutoff = min_abs_cutoff, mode = method,
                              return_details = return_details)
  # Return consistent list format
  if (is.list(result) && !is.data.frame(result)) {
    result
  } else {
    # When mode="zero", remove features that became all zeros
    if (method == "zero") {
      sample_cols <- as.matrix(result[, -1, drop = FALSE])
      keep_features <- rowSums(sample_cols) > 0
      result <- result[keep_features, , drop = FALSE]
    }
    list(table = result, detailed_leakage = NULL)
  }
}

#' Apply sparsity elbow-based filtering#'#' Detects elbow point in richness-depth curve and uses it to set#' a coverage cutoff for sample filtering.#'#' @param table Feature table (data.frame)#' @param detect Whether to run detection#' @param method Elbow detection method#' @param apply Whether to apply the filter#' @param multiplier MAD multiplier for applying elbow cutoff#' @param verbose Logical. Print progress messages?#'#' @return A list with filtered_table and elbow_result
#' @noRd
.apply_sparsity_elbow_filter <- function(table, detect, method, apply, multiplier, verbose = TRUE) {
  result <- list(filtered_table = table, elbow_result = NULL)

  if (!detect) {
    if (verbose) cat("Skipping sparsity elbow detection\n")
    return(result)
  }

  if (verbose) cat("Running sparsity elbow detection...\n")

  # Run elbow detection
  elbow_result <- identify_sparsity_elbow(table, method = method)
  result$elbow_result <- elbow_result

  if (!apply) {
    if (verbose) cat("Sparsity elbow detection completed but not applied\n")
    return(result)
  }

  # Use elbow threshold as coverage cutoff
  elbow_threshold <- elbow_result$kneedle_threshold
  if (is.null(elbow_threshold) || is.na(elbow_threshold)) {
    warning("Could not determine elbow threshold, skipping application")
    return(result)
  }

  # Apply MAD-based adjustment
  adjusted_cutoff <- elbow_threshold - multiplier * mad(elbow_result$richness_values)
  adjusted_cutoff <- max(adjusted_cutoff, 0)

  if (verbose) cat(sprintf("Applying sparsity elbow filter (cutoff: %.0f reads)...\n", adjusted_cutoff))
  result$filtered_table <- filter_by_coverage(table, min_reads = adjusted_cutoff)
  result
}

#' Apply depth-sparsity outlier filtering#'#' Identifies and optionally removes samples that are outliers in the#' depth-sparsity relationship.#'#' @param table Feature table (data.frame)#' @param detect Whether to run analysis#' @param metric Metric to analyze: "sparsity" or "richness"#' @param method Outlier detection method: "mad", "iqr", or "both"#' @param multiplier MAD/IQR multiplier#' @param direction Direction to flag: "high_sparsity", "low_sparsity", or "both"#' @param apply Whether to apply the filter#' @param verbose Logical. Print progress messages?#'#' @return A list with filtered_table and analysis_result
#' @noRd
.apply_depth_sparsity_filter <- function(table, detect, metric, method, multiplier,
                                          direction, apply, verbose = TRUE) {
  result <- list(filtered_table = table, analysis_result = NULL)

  if (!detect) {
    if (verbose) cat("Skipping depth-sparsity outlier analysis\n")
    return(result)
  }

  if (verbose) cat(sprintf("Running depth-sparsity analysis (metric: %s, method: %s)...\n",
                           metric, method))

  # Run analysis
  analysis_result <- analyze_depth_sparsity(table, metric = metric, outlier_method = method,
                                             multiplier = multiplier, direction = direction)
  result$analysis_result <- analysis_result

  if (!apply) {
    if (verbose) cat("Depth-sparsity analysis completed but not applied\n")
    return(result)
  }

  # Remove outliers
  if (verbose) cat("Removing depth-sparsity outliers...\n")
  result$filtered_table <- filter_depth_sparsity_outliers(table, metric = metric,
                                                           outlier_method = method,
                                                           multiplier = multiplier,
                                                           direction = direction)
  result
}

#' Apply abundance filtering#'#' Wrapper that handles different abundance filtering methods.#'#' @param table Feature table (data.frame)#' @param method Method: "none", "absolute", "relative", "relative_cutoff", "joint"#' @param threshold Threshold value#' @param min_samples Minimum number of samples#' @param logic Logical operator for joint filtering: "OR" or "AND"#' @param prevalence_threshold Prevalence threshold for joint filtering#' @param min_coverage_for_relative Minimum coverage for relative_cutoff method#' @param remove_features Whether to remove features (vs. setting to zero)#' @param verbose Logical. Print progress messages?#'#' @return Filtered feature table
#' @noRd
.apply_abundance_filter <- function(table, method, threshold, min_samples,
                                     logic, prevalence_threshold,
                                     min_coverage_for_relative, remove_features, verbose = TRUE) {
  if (method == "none") {
    if (verbose) cat("Skipping abundance filtering (method = 'none')\n")
    return(table)
  }

  if (verbose) cat(sprintf("Applying abundance filter (method: %s, threshold: %s)...\n",
                           method, threshold))

  switch(method,
    "absolute" = filter_features_by_abundance(
      table, threshold = threshold, mode = "absolute",
      min_samples = min_samples, remove_zeros = !remove_features
    ),
    "relative" = filter_features_by_abundance(
      table, threshold = threshold, mode = "relative",
      min_samples = min_samples, remove_zeros = !remove_features
    ),
    "relative_cutoff" = filter_by_relative_cutoff(
      table, min_coverage = min_coverage_for_relative,
      relative_threshold = threshold, remove_features = remove_features
    ),
    "joint" = filter_features_joint(
      table, abundance_threshold = threshold,
      prevalence_threshold = prevalence_threshold,
      mode = if (threshold > 1) "absolute" else "relative",
      logic = logic, remove_zeros = !remove_features
    )$table,
    stop("Unknown abundance method: ", method)
  )
}
