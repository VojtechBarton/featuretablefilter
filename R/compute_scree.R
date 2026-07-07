#' Compute scree/saturation diagnostics for filtering thresholds
#'
#' Evaluates feature and sample retention across a gradient of thresholds to help
#' users visualize how quickly their data collapses and choose an optimal cutoff.
#' This is particularly useful for MAD/IQR multiplier sweeps, coverage target sweeps
#' (Good's or Chao's), absolute count thresholds, or relative abundance cutoffs.
#'
#' The function generates diagnostic curves showing:
#' - Feature retention rate (% features kept at each threshold)
#' - Sample retention rate (% samples kept at each threshold)
#' - Total read retention (% reads preserved)
#' - Sparsity changes (proportion of zeros in the table)
#' - "Collapse rate" - derivative showing where rapid losses occur
#'
#' @param table A feature table (data.frame or matrix) with features as rows and samples as columns.
#'              First column should be feature IDs, remaining columns are sample counts.
#' @param type Type of threshold sweep to perform:
#'             \itemize{
#'               \item \code{"mad_multiplier"} - Sweep MAD multipliers (e.g., 1 to 5 for coverage filtering)
#'               \item \code{"iqr_multiplier"} - Sweep IQR multipliers (e.g., 0.5 to 3 for Tukey fence-based filtering)
#'               \item \code{"good_coverage"} - Sweep Good's coverage targets (0.70 to 0.99)
#'               \item \code{"chao_coverage"} - Sweep Chao's coverage targets (0.70 to 0.99)
#'               \item \code{"singleton_ratio"} - Sweep singleton ratio thresholds (0.01 to 0.30 for PCR artifact detection)
#'               \item \code{"cross_talk"} - Sweep cross-talk relative thresholds (0.0001 to 0.01 for index hopping)
#'               \item \code{"absolute_feature"} - Sweep absolute feature abundance thresholds
#'               \item \code{"relative_feature"} - Sweep relative abundance thresholds (0.00001 to 0.01)
#'               \item \code{"custom"} - Use custom threshold values provided in \code{thresholds}
#'             }
#' @param thresholds For \code{type = "custom"}, a numeric vector of threshold values to evaluate.
#'                   Ignored for other types (generated automatically).
#' @param n_steps Number of threshold steps to evaluate. Default is 20. More steps give smoother curves.
#' @param min_samples Minimum number of samples required for feature retention (passed to
#'                    \code{\link{filter_features_by_abundance}}). Default is 1.
#' @param verbose Logical. Print progress information. Default is TRUE.
#'
#' @return A list containing:
#'   \item{results}{data.frame with columns: threshold, n_features_retained, pct_features_retained,
#'                  n_samples_retained, pct_samples_retained, n_reads_retained, pct_reads_retained,
#'                  sparsity, collapse_rate}
#'   \item{summary}{list with key statistics including elbow points and saturation metrics}
#'   \item{type}{the type of sweep performed}
#'   \item{parameters}{list of parameters used}
#'
#' @export
#'
#' @examples
#' # Sweep MAD multipliers from 1 to 5 (for coverage-based sample filtering)
#' # scree <- compute_scree(my_table, type = "mad_multiplier", n_steps = 20)
#' # plot(scree$results$threshold, scree$results$pct_features_retained, type = "l")
#'
#' # Sweep IQR multipliers from 0.5 to 3 (Tukey fence-based sample filtering)
#' # scree <- compute_scree(my_table, type = "iqr_multiplier", n_steps = 20)
#'
#' # Sweep Good's coverage targets from 70% to 99%
#' # scree <- compute_scree(my_table, type = "good_coverage", n_steps = 20)
#'
#' # Sweep Chao's coverage targets from 70% to 99%
#' # scree <- compute_scree(my_table, type = "chao_coverage", n_steps = 20)
#'
#' # Sweep singleton ratio thresholds from 1% to 30% (for PCR artifact detection)
#' # scree <- compute_scree(my_table, type = "singleton_ratio", n_steps = 20)
#'
#' # Sweep cross-talk relative thresholds from 0.01% to 1% (for index hopping detection)
#' # scree <- compute_scree(my_table, type = "cross_talk", n_steps = 20)
#'
#' # Sweep absolute feature abundance thresholds
#' # scree <- compute_scree(my_table, type = "absolute_feature", n_steps = 30)
#'
#' # Custom threshold sweep
#' # scree <- compute_scree(my_table, type = "custom", thresholds = c(1, 2, 5, 10, 20, 50))
compute_scree <- function(table, type = c("mad_multiplier", "iqr_multiplier",
                                           "good_coverage", "chao_coverage",
                                           "singleton_ratio", "cross_talk",
                                           "absolute_feature", "relative_feature", "custom"),
                          thresholds = NULL, n_steps = 20, min_samples = 1,
                          count_type = "both", verbose = TRUE) {
  # Validate inputs
  if (!is.data.frame(table) && !is.matrix(table)) {
    stop("table must be a data.frame or matrix")
  }

  if (ncol(table) < 2) {
    stop("table must have at least one sample column")
  }

  type <- match.arg(type)

  if (type == "custom" && is.null(thresholds)) {
    stop("thresholds must be provided when type = 'custom'")
  }

  # Extract abundance data
  abundances <- as.matrix(table[, -1, drop = FALSE])
  feature_ids <- table[, 1, drop = FALSE]

  # Calculate baseline metrics
  total_reads_baseline <- sum(abundances)
  n_features_baseline <- nrow(abundances)
  n_samples_baseline <- ncol(abundances)
  n_zeros_baseline <- sum(abundances == 0)
  baseline_sparsity <- n_zeros_baseline / (n_features_baseline * n_samples_baseline)

  # Generate threshold sequence based on type
  if (type == "mad_multiplier") {
    # Sweep MAD multipliers for coverage-based filtering
    thresholds <- seq(1, 5, length.out = n_steps)
  } else if (type == "iqr_multiplier") {
    # Sweep IQR multipliers for Tukey fence-based filtering
    thresholds <- seq(0.5, 3, length.out = n_steps)
  } else if (type == "good_coverage") {
    # Sweep Good's coverage targets (70% to 99%)
    thresholds <- seq(0.70, 0.99, length.out = n_steps)
  } else if (type == "chao_coverage") {
    # Sweep Chao's coverage targets (70% to 99%)
    thresholds <- seq(0.70, 0.99, length.out = n_steps)
  } else if (type == "singleton_ratio") {
    # Sweep singleton ratio thresholds (1% to 30% for PCR artifact detection)
    thresholds <- seq(0.01, 0.30, length.out = n_steps)
  } else if (type == "cross_talk") {
    # Sweep cross-talk relative thresholds (0.01% to 1% of feature max)
    thresholds <- seq(0.0001, 0.01, length.out = n_steps)
  } else if (type == "absolute_feature") {
    # Sweep absolute abundance (log-spaced for better resolution at low values)
    # Only generate thresholds if not provided by user
    if (is.null(thresholds)) {
      max_abund <- max(abundances)
      if (max_abund < n_steps) {
        thresholds <- 1:max(min(max_abund, 50), ceiling(n_steps/2))
      } else {
        thresholds <- unique(c(1:ceiling(n_steps/2),
                               round(exp(seq(log(2), log(max_abund), length.out = n_steps/2)))))
      }
      thresholds <- sort(unique(thresholds))
    }
  } else if (type == "relative_feature") {
    # Sweep relative abundance (0.001% to 1%)
    thresholds <- seq(0.00001, 0.01, length.out = n_steps)
  }
  # For custom, use provided thresholds (defaults to absolute feature filtering behavior)
  # Note: custom type applies thresholds as absolute feature abundance cutoffs

  # Initialize filtered_abundances to avoid "object not found" error
  filtered_abundances <- abundances

  if (verbose) {
    message(sprintf("Evaluating %d threshold values...", length(thresholds)))
  }

  # Initialize results matrix
  results <- matrix(NA, nrow = length(thresholds), ncol = 9)
  colnames(results) <- c("threshold", "n_features_retained", "pct_features_retained",
                         "n_samples_retained", "pct_samples_retained",
                         "n_reads_retained", "pct_reads_retained",
                         "sparsity", "collapse_rate")

  # Evaluate each threshold
  for (i in seq_along(thresholds)) {
    thresh <- thresholds[i]

    if (type == "mad_multiplier") {
      # Apply coverage-based sample filtering using MAD multiplier
      sample_sums <- colSums(abundances)
      median_cov <- median(sample_sums)
      mad_cov <- mad(sample_sums)
      cutoff <- median_cov - thresh * mad_cov
      cutoff <- max(cutoff, 0)

      keep_samples <- sample_sums >= cutoff
      filtered_abundances <- abundances[, keep_samples, drop = FALSE]

      # All features retained (we're filtering samples, not features)
      n_features_retained <- n_features_baseline
      n_samples_retained <- sum(keep_samples)
      n_reads_retained <- sum(filtered_abundances)

    } else if (type == "iqr_multiplier") {
      # Apply coverage-based sample filtering using IQR multiplier (Tukey fences)
      sample_sums <- colSums(abundances)
      q1 <- quantile(sample_sums, 0.25)
      q3 <- quantile(sample_sums, 0.75)
      iqr <- q3 - q1
      cutoff <- q1 - thresh * iqr
      cutoff <- max(cutoff, 0)

      keep_samples <- sample_sums >= cutoff
      filtered_abundances <- abundances[, keep_samples, drop = FALSE]

      n_features_retained <- n_features_baseline
      n_samples_retained <- sum(keep_samples)
      n_reads_retained <- sum(filtered_abundances)

    } else if (type == "good_coverage") {
      # Apply Good's coverage-based sample filtering
      sample_totals <- colSums(abundances)
      singletons_per_sample <- apply(abundances, 2, function(x) sum(x == 1))
      sample_coverage <- 1 - (singletons_per_sample / sample_totals)
      sample_coverage[is.na(sample_coverage)] <- 0

      keep_samples <- sample_coverage >= thresh
      filtered_abundances <- abundances[, keep_samples, drop = FALSE]

      n_features_retained <- n_features_baseline
      n_samples_retained <- sum(keep_samples)
      n_reads_retained <- sum(filtered_abundances)

    } else if (type == "chao_coverage") {
      # Apply Chao's coverage-based sample filtering
      sample_totals <- colSums(abundances)
      n_features_per_sample <- colSums(abundances > 0)
      singletons_per_sample <- apply(abundances, 2, function(x) sum(x == 1))
      doubletons_per_sample <- apply(abundances, 2, function(x) sum(x == 2))

      # Calculate Chao's coverage per sample using standard estimator formula:
      # C_hat = 1 - (f1/n) * [(n-1)*f1 / ((n-1)*f1 + 2*f2)]
      sample_coverage <- sapply(seq_along(sample_totals), function(i) {
        n <- sample_totals[i]
        f1 <- singletons_per_sample[i]
        f2 <- doubletons_per_sample[i]
        if (n == 0) return(0)
        denom <- (n - 1) * f1 + 2 * f2
        coverage <- if (denom == 0) 1 - f1 / n else 1 - (f1 / n) * (((n - 1) * f1) / denom)
        max(0, min(1, coverage))
      })

      keep_samples <- sample_coverage >= thresh
      filtered_abundances <- abundances[, keep_samples, drop = FALSE]

      n_features_retained <- n_features_baseline
      n_samples_retained <- sum(keep_samples)
      n_reads_retained <- sum(filtered_abundances)

    } else if (type == "singleton_ratio") {
      # Apply singleton ratio-based sample filtering
      sample_totals <- colSums(abundances)
      singleton_counts <- apply(abundances, 2, function(x) sum(x == 1))
      doubleton_counts <- apply(abundances, 2, function(x) sum(x == 2))

      if (count_type == "singleton") {
        low_count_sum <- singleton_counts
      } else if (count_type == "doubleton") {
        low_count_sum <- doubleton_counts
      } else {  # "both"
        low_count_sum <- singleton_counts + doubleton_counts
      }

      # Calculate ratio of low-count reads to total reads per sample
      ratio_vector <- ifelse(sample_totals > 0, low_count_sum / sample_totals, NA)

      # Keep samples with ratio below threshold (or NA due to zero reads)
      keep_samples <- is.na(ratio_vector) | (ratio_vector <= thresh)
      filtered_abundances <- abundances[, keep_samples, drop = FALSE]

      n_features_retained <- n_features_baseline
      n_samples_retained <- sum(keep_samples)
      n_reads_retained <- sum(filtered_abundances)

    } else if (type == "absolute_feature") {
      # Apply absolute feature abundance filtering
      keep_features <- rowSums(abundances >= thresh) >= min_samples
      filtered_abundances <- abundances[keep_features, , drop = FALSE]

      n_features_retained <- sum(keep_features)
      n_samples_retained <- n_samples_baseline  # No sample filtering
      n_reads_retained <- sum(filtered_abundances)

    } else if (type == "cross_talk") {
      # Apply cross-talk/index-hopping filtering
      # For scree analysis, we count how many features would be removed
      # when using "remove_feature" mode at each threshold
      feature_max <- apply(abundances, 1, max)
      dynamic_threshold <- outer(feature_max, rep(1, ncol(abundances))) * thresh

      # Identify leakage: values > 0 AND < dynamic_threshold
      leakage_mask <- abundances > 0 & abundances < dynamic_threshold

      # Count features with any leakage (these would be removed in remove_feature mode)
      features_with_leakage <- rowSums(leakage_mask) > 0
      keep_features <- !features_with_leakage

      filtered_abundances <- abundances[keep_features, , drop = FALSE]

      n_features_retained <- sum(keep_features)
      n_samples_retained <- n_samples_baseline
      n_reads_retained <- sum(filtered_abundances)

    } else if (type == "relative_feature") {
      # Apply relative abundance filtering
      sample_totals <- colSums(abundances)
      rel_abundances <- sweep(abundances, 2, sample_totals, FUN = "/")
      keep_features <- rowSums(rel_abundances >= thresh) >= min_samples
      filtered_abundances <- abundances[keep_features, , drop = FALSE]

      n_features_retained <- sum(keep_features)
      n_samples_retained <- n_samples_baseline
      n_reads_retained <- sum(filtered_abundances)
    } else if (type == "custom") {
      # Custom thresholds - treat as absolute feature abundance filtering
      keep_features <- rowSums(abundances >= thresh) >= min_samples
      filtered_abundances <- abundances[keep_features, , drop = FALSE]

      n_features_retained <- sum(keep_features)
      n_samples_retained <- n_samples_baseline
      n_reads_retained <- sum(filtered_abundances)
    }

    # Handle edge case of empty result
    if (nrow(filtered_abundances) == 0 || ncol(filtered_abundances) == 0) {
      n_features_retained <- 0
      n_samples_retained <- 0
      n_reads_retained <- 0
      sparsity <- 1
    } else {
      n_zeros_filtered <- sum(filtered_abundances == 0)
      sparsity <- n_zeros_filtered / (nrow(filtered_abundances) * ncol(filtered_abundances))
    }

    # Store results
    results[i, ] <- c(
      thresh,
      n_features_retained,
      (n_features_retained / n_features_baseline) * 100,
      n_samples_retained,
      (n_samples_retained / n_samples_baseline) * 100,
      n_reads_retained,
      (n_reads_retained / total_reads_baseline) * 100,
      sparsity,
      NA  # collapse_rate calculated after loop
    )

    if (verbose && i %% max(1, length(thresholds) %/% 5) == 0) {
      message(sprintf("  Completed %d/%d (threshold = %.4f)",
                      i, length(thresholds), thresh))
    }
  }

  # Calculate collapse rate (negative derivative of retention)
  # Higher values indicate faster collapse at that threshold
  retention_pct <- results[, "pct_features_retained"]
  if (length(retention_pct) > 1) {
    delta_retention <- diff(retention_pct)
    delta_thresh <- diff(as.numeric(results[, "threshold"]))
    collapse_rate <- -delta_retention / abs(delta_thresh + .Machine$double.eps)
    # Pad with NA to match length
    collapse_rate <- c(NA, collapse_rate)
  } else {
    collapse_rate <- NA
  }
  results[, "collapse_rate"] <- collapse_rate

  # Convert to data.frame
  results_df <- as.data.frame(results)

  # Identify elbow points (where collapse rate peaks)
  valid_rates <- results_df$collapse_rate[!is.na(results_df$collapse_rate)]
  if (length(valid_rates) > 0) {
    max_rate_idx <- which.max(results_df$collapse_rate[!is.na(results_df$collapse_rate)]) + 1
    elbow_threshold <- results_df$threshold[max_rate_idx]
    elbow_retention <- results_df$pct_features_retained[max_rate_idx]
  } else {
    elbow_threshold <- NA
    elbow_retention <- NA
  }

  # Calculate saturation metrics
  # How much does the curve flatten out at high thresholds?
  if (n_steps >= 3) {
    last_three_retention <- tail(results_df$pct_features_retained, 3)
    saturation_slope <- mean(diff(last_three_retention))
    saturated <- abs(saturation_slope) < 1  # Less than 1% change per step
  } else {
    saturation_slope <- NA
    saturated <- NA
  }

  summary_list <- list(
    n_thresholds_evaluated = length(thresholds),
    threshold_range = range(thresholds),
    elbow_point = list(
      threshold = elbow_threshold,
      retention_at_elbow = elbow_retention
    ),
    saturation = list(
      final_retention = results_df$pct_features_retained[length(results_df$pct_features_retained)],
      final_threshold = results_df$threshold[length(results_df$threshold)],
      saturation_slope = saturation_slope,
      is_saturated = saturated
    ),
    baseline = list(
      n_features = n_features_baseline,
      n_samples = n_samples_baseline,
      total_reads = total_reads_baseline,
      sparsity = baseline_sparsity
    )
  )

  if (verbose) {
    message(sprintf("\nScree analysis complete."))
    message(sprintf("  Elbow point: threshold = %.4f, retention = %.1f%%",
                    elbow_threshold, elbow_retention))
    message(sprintf("  Final retention: %.1f%% at threshold = %.4f",
                    summary_list$saturation$final_retention,
                    summary_list$saturation$final_threshold))
  }

  return(list(
    results = results_df,
    summary = summary_list,
    type = type,
    parameters = list(n_steps = n_steps, min_samples = min_samples)
  ))
}


#' Plot scree/saturation diagnostic results
#'
#' Creates a multi-panel visualization of scree analysis results showing
#' retention curves, sparsity changes, and collapse rates.
#'
#' @param scree_obj Output from \code{\link{compute_scree}}
#' @param main Main title for the plot
#' @param show_collapse Logical. Show collapse rate panel. Default is TRUE.
#' @param color Theme color palette. Default is "blue". Options: "blue", "green", "red", "purple"
#'
#' @return A grid object containing the plot (requires ggplot2)
#'
#' @export
#'
#' @examples
#' # scree <- compute_scree(my_table, type = "absolute_feature")
#' # plot_scree(scree)
plot_scree <- function(scree_obj, main = "Filtering Threshold Scree Analysis",
                       show_collapse = TRUE, color = "blue", verbose = TRUE, ...) {
  # Check for ggplot2
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting. Please install it.")
  }

  library(ggplot2)

  results <- scree_obj$results

  # Color mapping
  color_map <- list(
    blue = "#2E86AB",
    green = "#239B56",
    red = "#C0392B",
    purple = "#884EA0"
  )
  primary_color <- color_map[[color]]
  if (is.null(primary_color)) primary_color <- color_map[["blue"]]

  # Create base theme
  base_theme <- theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      panel.grid.minor = element_line(color = "grey90"),
      strip.background = element_rect(fill = "grey95")
    )

  # Panel 1: Retention curves (features, samples, reads)
  p1 <- ggplot(results, aes(x = threshold)) +
    geom_line(aes(y = pct_features_retained, color = "Features"), linewidth = 1) +
    geom_line(aes(y = pct_samples_retained, color = "Samples"), linewidth = 1) +
    geom_line(aes(y = pct_reads_retained, color = "Reads"), linewidth = 1) +
    scale_color_manual(values = c(
      "Features" = primary_color,
      "Samples" = "#E74C3C",
      "Reads" = "#27AE60"
    )) +
    labs(
      title = "Retention Curves",
      x = "Threshold",
      y = "Retention (%)",
      color = "Metric"
    ) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
    base_theme

  # Panel 2: Sparsity changes
  p2 <- ggplot(results, aes(x = threshold, y = sparsity * 100)) +
    geom_line(color = primary_color, linewidth = 1) +
    geom_area(alpha = 0.3, fill = primary_color) +
    labs(
      title = "Sparsity Change",
      x = "Threshold",
      y = "Sparsity (% zeros)"
    ) +
    base_theme

  # Panel 3: Collapse rate (derivative)
  if (show_collapse && any(!is.na(results$collapse_rate))) {
    p3 <- ggplot(results, aes(x = threshold, y = collapse_rate)) +
      geom_line(color = "#C0392B", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      labs(
        title = "Collapse Rate (Derivative)",
        x = "Threshold",
        y = "Rate of Loss (% per unit)"
      ) +
      base_theme

    panels <- list(p1, p2, p3)
  } else {
    panels <- list(p1, p2)
  }

  # Combine panels
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    combined <- wrap_plots(panels, ncol = 1, heights = c(1.2, 1, 0.8)) +
      plot_annotation(title = main, theme = theme(plot.title = element_text(size = 16)))
  } else {
    # Fallback without patchwork
    combined <- p1
    warning("patchwork recommended for multi-panel layout. Install with install.packages('patchwork').")
  }

  return(combined)
}
