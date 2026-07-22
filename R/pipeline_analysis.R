# =============================================================================
# Pipeline Analysis Step Wrapper Functions (Internal)
# =============================================================================
# These functions wrap the existing analysis functions for use in the pipeline.
# All functions are internal (. prefix).
# =============================================================================

#' Run sparsity elbow detection (diagnostic only)#'#' Wrapper around identify_sparsity_elbow() for diagnostic purposes.#'#' @param table Feature table (data.frame)#' @param method Elbow detection method: "kneedle", "max_derivative", or "second_derivative"#' @param verbose Logical. Print progress messages?#'#' @return Result from identify_sparsity_elbow()
#' @noRd
.run_sparsity_elbow_analysis <- function(table, method, verbose = TRUE) {
  if (verbose) cat(sprintf("Running sparsity elbow detection (method: %s)...\n", method))
  identify_sparsity_elbow(table, method = method)
}

#' Run depth-sparsity outlier analysis (diagnostic only)#'#' Wrapper around analyze_depth_sparsity() for diagnostic purposes.#'#' @param table Feature table (data.frame)#' @param metric Metric to analyze: "sparsity" or "richness"#' @param method Outlier detection method: "mad", "iqr", or "both"#' @param multiplier MAD/IQR multiplier#' @param direction Direction to flag: "high_sparsity", "low_sparsity", or "both"#' @param verbose Logical. Print progress messages?#'#' @return Result from analyze_depth_sparsity()
#' @noRd
.run_depth_sparsity_analysis <- function(table, metric, method, multiplier, direction, verbose = TRUE) {
  if (verbose) cat(sprintf("Running depth-sparsity analysis (metric: %s, method: %s)...\n",
                           metric, method))
  analyze_depth_sparsity(table, metric = metric, outlier_method = method,
                          multiplier = multiplier, direction = direction)
}

#' Run scree/saturation analysis#'#' Wrapper around compute_scree() for comprehensive filtering threshold exploration.#'#' @param table Original feature table (data.frame)#' @param type Scree sweep type: "mad_multiplier", "absolute_feature", "relative_feature", or "custom"#' @param n_steps Number of steps for the sweep#' @param thresholds Custom thresholds vector (for type = "custom")#' @param verbose Logical. Print progress messages?#'#' @return Result from compute_scree()
#' @noRd
.run_scree_analysis <- function(table, type, n_steps, thresholds, verbose = TRUE) {
  if (verbose) cat(sprintf("Running scree analysis (type: %s, steps: %d)...\n", type, n_steps))

  switch(type,
    "mad_multiplier" = compute_scree(
      table, type = "mad_multiplier",
      n_steps = n_steps
    ),
    "absolute_feature" = compute_scree(
      table, type = "absolute_feature",
      thresholds = seq(1, n_steps),
      n_steps = n_steps
    ),
    "relative_feature" = compute_scree(
      table, type = "relative_feature",
      thresholds = seq(0.1, 1, length.out = n_steps),
      n_steps = n_steps
    ),
    "custom" = compute_scree(
      table, type = "custom",
      thresholds = thresholds
    ),
    stop("Unknown scree type: ", type)
  )
}

#' Compute table statistics#'#' Helper function to get basic statistics about a feature table.#'#' @param table Feature table (data.frame)#'#' @return A list with features, samples, and reads counts
#' @noRd
.get_table_stats <- function(table) {
  n_features <- nrow(table)
  n_samples <- ncol(table) - 1  # Exclude feature_id column
  n_reads <- sum(as.matrix(table[, -1, drop = FALSE]))
  list(features = n_features, samples = n_samples, reads = n_reads)
}
