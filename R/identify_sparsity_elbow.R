#' Identify the sparsity elbow in sample coverage distribution
#'
#' Detects the critical sequencing depth threshold where low-prevalence ASVs
#' begin to crash, indicating insufficient sampling depth. This method analyzes
#' the structural stability of the feature table by examining how observed
#' richness changes with sequencing depth, rather than relying on singleton counts.
#'
#' The algorithm:
#' 1. Ranks samples by sequencing depth (total reads)
#' 2. Calculates cumulative feature discovery curve
#' 3. Computes the rate of new feature discovery at each depth level
#' 4. Identifies the "elbow" point where richness accumulation sharply declines
#'
#' Samples below this elbow threshold are likely undersampled and may contain
#' mostly noise/spurious features rather than true biological signal.
#'
#' @param table A feature table (data.frame or matrix) with features as rows and samples as columns.
#'              First column should be feature IDs, remaining columns are sample counts.
#' @param method Elbow detection method:
#'               \itemize{
#'                 \item \code{"kneedle"} - Kneedle algorithm for automatic elbow detection
#'                 \item{\code{"max_derivative"}} - Point of maximum negative derivative in richness curve
#'                 \item{\code{"second_derivative"}} - Point where second derivative peaks (inflection point)
#'               }
#' @param smooth_window Window size for smoothing the richness curve before derivative calculation.
#'                      Default is 5. Set to 1 for no smoothing.
#' @param min_samples Minimum number of samples required for analysis. Default is 5.
#' @param percentile_range Range of percentiles to consider for elbow detection (lower, upper).
#'                         Default is c(10, 90) to avoid extreme outliers.
#'
#' @return A list containing:
#'   \item{elbow_threshold}{The identified sequencing depth threshold}
#'   \item{samples_above_elbow}{Number of samples above the threshold}
#'   \item{samples_below_elbow}{Number of samples below the threshold}
#'   \item{richness_curve}{data.frame with sample_rank, depth, observed_richness, smoothed_richness, derivative}
#'   \item{recommendation}{Character string with filtering recommendation}
#'   \item{metrics}{List of diagnostic metrics including R-squared of fit, curvature at elbow}
#'
#' @export
#'
#' @examples
#' # Identify sparsity elbow in your data
#' # elbow_result <- identify_sparsity_elbow(my_table)
#' # filter_cutoff <- elbow_result$elbow_threshold
#' # filtered_samples <- filter_by_coverage(my_table, min_reads = filter_cutoff)
#'
#' # Use different detection methods
#' # elbow_kneedle <- identify_sparsity_elbow(my_table, method = "kneedle")
#' # elbow_derivative <- identify_sparsity_elbow(my_table, method = "max_derivative")
identify_sparsity_elbow <- function(table, method = c("kneedle", "max_derivative", "second_derivative"),
                                     smooth_window = 5, min_samples = 5,
                                     percentile_range = c(10, 90)) {
  # Validate inputs
  if (!is.data.frame(table) && !is.matrix(table)) {
    stop("table must be a data.frame or matrix")
  }

  if (ncol(table) < 2) {
    stop("table must have at least one sample column")
  }

  method <- match.arg(method)

  if (length(percentile_range) != 2 ||
      percentile_range[1] >= percentile_range[2] ||
      percentile_range[1] < 0 || percentile_range[2] > 100) {
    stop("percentile_range must be a vector of two values between 0 and 100, with first < second")
  }

  # Extract sample coverages
  sample_data <- table[, -1, drop = FALSE]
  sample_depths <- colSums(sample_data)
  sample_names <- colnames(sample_data)

  n_samples <- length(sample_depths)
  if (n_samples < min_samples) {
    stop(sprintf("Need at least %d samples for sparsity elbow analysis", min_samples))
  }

  # Calculate observed richness (number of non-zero features) per sample
  sample_richness <- colSums(sample_data > 0)

  # Create sorted dataframe by depth
  depth_order <- order(sample_depths, decreasing = TRUE)
  ranked_data <- data.frame(
    sample_name = sample_names[depth_order],
    depth = sample_depths[depth_order],
    richness = sample_richness[depth_order],
    rank = seq_along(sample_depths)
  )

  # Apply percentile range to focus on meaningful part of distribution
  lower_pct <- quantile(ranked_data$depth, probs = percentile_range[1] / 100)
  upper_pct <- quantile(ranked_data$depth, probs = percentile_range[2] / 100)
  valid_mask <- ranked_data$depth >= lower_pct & ranked_data$depth <= upper_pct

  # If too few points in range, expand to full range
  if (sum(valid_mask) < 3) {
    valid_mask <- rep(TRUE, n_samples)
  }

  valid_data <- ranked_data[valid_mask, ]

  # Smooth the richness curve
  if (smooth_window > 1 && nrow(valid_data) > smooth_window) {
    smoothed <- zoo::rollapply(valid_data$richness, width = smooth_window,
                                FUN = mean, align = "center", fill = NA, partial = TRUE)
    valid_data$smoothed_richness <- smoothed
  } else {
    valid_data$smoothed_richness <- valid_data$richness
  }

  # Normalize for elbow detection (scale to 0-1 range)
  depth_norm <- (valid_data$depth - min(valid_data$depth)) /
                (max(valid_data$depth) - min(valid_data$depth) + .Machine$double.eps)
  richness_norm <- (valid_data$smoothed_richness - min(valid_data$smoothed_richness)) /
                   (max(valid_data$smoothed_richness) - min(valid_data$smoothed_richness) + .Machine$double.eps)

  # Calculate derivatives (rate of change)
  # Negative derivative because we're going from high to low depth
  valid_data$derivative <- c(NA, diff(richness_norm) / (diff(depth_norm) + .Machine$double.eps))

  # Second derivative for inflection point detection
  valid_data$second_derivative <- c(NA, NA,
    diff(diff(richness_norm) / (diff(depth_norm) + .Machine$double.eps)) /
    (diff(depth_norm)[-1] + .Machine$double.eps))

  # Identify elbow based on method
  elbow_idx <- NA

  if (method == "kneedle") {
    # Kneedle algorithm: find point furthest from the line connecting endpoints
    # This is the point of maximum perpendicular distance from the chord
    x <- depth_norm
    y <- richness_norm

    # Line from (0, y_at_0) to (1, y_at_1)
    y_start <- y[which.min(x)]
    y_end <- y[which.max(x)]

    # Equation of line: y = y_start + (y_end - y_start) * x
    # Distance from point (x_i, y_i) to line: |y_i - (y_start + (y_end-y_start)*x_i)| / sqrt(1 + m^2)
    m <- (y_end - y_start) / (max(x) - min(x) + .Machine$double.eps)
    predicted_y <- y_start + m * x
    distances <- abs(y - predicted_y) / sqrt(1 + m^2)

    elbow_idx <- which.max(distances)

  } else if (method == "max_derivative") {
    # Find point with steepest negative slope (maximum magnitude of negative derivative)
    valid_deriv <- valid_data$derivative[!is.na(valid_data$derivative)]
    deriv_idx <- which(!is.na(valid_data$derivative))

    # Maximum negative derivative indicates where richness crashes
    elbow_local_idx <- deriv_idx[which.min(valid_deriv)]  # Most negative
    elbow_idx <- elbow_local_idx

  } else if (method == "second_derivative") {
    # Find inflection point where second derivative peaks
    valid_second <- valid_data$second_derivative[!is.na(valid_data$second_derivative)]
    second_idx <- which(!is.na(valid_data$second_derivative))

    if (length(valid_second) > 0) {
      elbow_local_idx <- second_idx[which.max(valid_second)]
      elbow_idx <- elbow_local_idx
    } else {
      # Fallback to max_derivative
      valid_deriv <- valid_data$derivative[!is.na(valid_data$derivative)]
      deriv_idx <- which(!is.na(valid_data$derivative))
      elbow_idx <- deriv_idx[which.min(valid_deriv)]
    }
  }

  # Extract elbow threshold
  if (!is.na(elbow_idx) && elbow_idx <= nrow(valid_data)) {
    elbow_threshold <- valid_data$depth[elbow_idx]
    elbow_richness <- valid_data$smoothed_richness[elbow_idx]
    elbow_depth_norm <- depth_norm[elbow_idx]
    elbow_curvature <- valid_data$derivative[elbow_idx]
  } else {
    # Fallback: use median depth
    elbow_threshold <- median(sample_depths)
    elbow_richness <- NA
    elbow_depth_norm <- NA
    elbow_curvature <- NA
  }

  # Count samples above/below elbow
  samples_above <- sum(sample_depths >= elbow_threshold)
  samples_below <- sum(sample_depths < elbow_threshold)

  # Calculate diagnostic metrics
  # R-squared of linear fit to the "stable" region (above elbow)
  stable_data <- valid_data[valid_data$depth >= elbow_threshold, ]
  if (nrow(stable_data) >= 2) {
    lm_fit <- lm(smoothed_richness ~ depth, data = stable_data)
    r_squared <- summary(lm_fit)$r.squared
  } else {
    r_squared <- NA
  }

  # Curvature ratio: how much steeper is the curve below vs above elbow
  above_deriv <- valid_data$derivative[valid_data$depth >= elbow_threshold & !is.na(valid_data$derivative)]
  below_deriv <- valid_data$derivative[valid_data$depth < elbow_threshold & !is.na(valid_data$derivative)]

  if (length(above_deriv) > 0 && length(below_deriv) > 0) {
    curvature_ratio <- abs(mean(below_deriv)) / (abs(mean(above_deriv)) + .Machine$double.eps)
  } else {
    curvature_ratio <- NA
  }

  # Generate recommendation
  if (samples_below / n_samples > 0.2) {
    recommendation <- sprintf(
      "Consider filtering out %d samples (%.1f%%) below depth %d reads. These samples show sharp richness decline indicating undersampling.",
      samples_below, (samples_below / n_samples) * 100, round(elbow_threshold)
    )
  } else if (curvature_ratio > 3) {
    recommendation <- sprintf(
      "Elbow detected at %d reads with %.1fx steeper decline below threshold. %d samples fall below this critical depth.",
      round(elbow_threshold), curvature_ratio, samples_below
    )
  } else {
    recommendation <- sprintf(
      "No sharp elbow detected. Richness-depth relationship is relatively stable. Current median depth (%d) may be adequate.",
      median(sample_depths)
    )
  }

  # Build richness curve for full dataset (not just valid range)
  ranked_data$smoothed_richness <- NA
  ranked_data$derivative <- NA
  ranked_data$second_derivative <- NA

  # Map back smoothed values
  for (i in seq_len(nrow(valid_data))) {
    orig_idx <- which(ranked_data$sample_name == valid_data$sample_name[i])
    if (length(orig_idx) > 0) {
      ranked_data$smoothed_richness[orig_idx] <- valid_data$smoothed_richness[i]
      ranked_data$derivative[orig_idx] <- valid_data$derivative[i]
      ranked_data$second_derivative[orig_idx] <- valid_data$second_derivative[i]
    }
  }

  return(list(
    elbow_threshold = elbow_threshold,
    samples_above_elbow = samples_above,
    samples_below_elbow = samples_below,
    richness_curve = ranked_data,
    recommendation = recommendation,
    metrics = list(
      r_squared_stable_region = r_squared,
      curvature_ratio = curvature_ratio,
      elbow_depth_normalized = elbow_depth_norm,
      elbow_curvature = elbow_curvature,
      method_used = method,
      percentile_range = percentile_range,
      smooth_window = smooth_window
    )
  ))
}


#' Plot sparsity elbow analysis results
#'
#' Creates visualization of the richness-depth relationship with elbow point marked.
#'
#' @param elbow_result Output from \code{\link{identify_sparsity_elbow}}
#' @param main Main title for the plot
#'
#' @return A ggplot object
#'
#' @export
#'
#' @examples
#' # elbow <- identify_sparsity_elbow(my_table)
#' # plot_sparsity_elbow(elbow)
plot_sparsity_elbow <- function(elbow_result, main = "Sparsity Elbow Analysis") {
  # Check for ggplot2
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting. Please install it.")
  }

  curve_df <- elbow_result$richness_curve
  elbow_thresh <- elbow_result$elbow_threshold

  base_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = ggplot2::element_text(size = 11),
      axis.text = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_line(color = "grey90")
    )

  # Panel 1: Richness vs Depth with elbow marked
  p1 <- ggplot2::ggplot(curve_df, ggplot2::aes(x = depth, y = richness)) +
    ggplot2::geom_point(alpha = 0.5, color = "grey50", size = 1.5) +
    ggplot2::geom_line(ggplot2::aes(y = smoothed_richness), color = "#2E86AB", linewidth = 1) +
    ggplot2::geom_vline(xintercept = elbow_thresh, color = "#C0392B", linetype = "dashed", linewidth = 1) +
    ggplot2::geom_hline(yintercept = NA, color = "#C0392B", linetype = "dotted", alpha = 0.5) +
    ggplot2::annotate("text", x = elbow_thresh * 1.05, y = max(curve_df$richness) * 0.95,
             label = paste0("Elbow: ", round(elbow_thresh), " reads"),
             color = "#C0392B", hjust = 0, size = 4) +
    ggplot2::labs(
      title = "Sample Richness vs Sequencing Depth",
      subtitle = paste0("Samples above elbow: ", elbow_result$samples_above_elbow,
                        " | Below elbow: ", elbow_result$samples_below_elbow),
      x = "Sequencing Depth (reads)",
      y = "Observed Richness (ASVs)"
    ) +
    scale_x_log10() +
    scale_y_log10() +
    base_theme

  # Panel 2: Derivative (rate of change) showing where crash occurs
  deriv_df <- curve_df[!is.na(curve_df$derivative), ]
  if (nrow(deriv_df) > 0) {
    p2 <- ggplot(deriv_df, aes(x = depth, y = derivative)) +
      geom_line(color = "#884EA0", linewidth = 1) +
      geom_vline(xintercept = elbow_thresh, color = "#C0392B", linetype = "dashed", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.5) +
      labs(
        title = "Rate of Richness Change (Derivative)",
        subtitle = "Negative spikes indicate rapid richness loss at lower depths",
        x = "Sequencing Depth (reads)",
        y = "Normalized Richness Derivative"
      ) +
      scale_x_log10() +
      base_theme
  } else {
    p2 <- ggplot() +
      geom_text(x = 0.5, y = 0.5, label = "Insufficient data for derivative plot") +
      theme_void()
  }

  # Combine panels
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- (p1 / p2) +
      patchwork::plot_annotation(title = main,
                      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16)))
  } else {
    combined <- p1
    warning("patchwork recommended for multi-panel layout. Install with install.packages('patchwork').")
  }

  return(combined)
}
