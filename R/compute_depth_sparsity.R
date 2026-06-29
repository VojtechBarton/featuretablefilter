#' Compute depth-sparsity relationship and identify outlier samples
#'
#' Analyzes the relationship between sequencing depth (total reads) and sparsity
#' (proportion of zeros) or observed richness per sample. Uses MAD/IQR-based
# ' outlier detection to flag samples that deviate significantly from the expected
#' depth-sparsity curve, which may indicate poor quality or technical artifacts.
#'
#' The function:
#' 1. Calculates depth and sparsity/richness for each sample
#' 2. Fits a robust regression line to the depth-sparsity relationship
#' 3. Identifies outliers using MAD or IQR-based methods
#' 4. Provides recommendations for which samples to consider removing
#'
#' @param table A feature table (data.frame or matrix) with features as rows and samples as columns.
#'              First column should be feature IDs, remaining columns are sample counts.
#' @param metric What to plot against depth: \code{"sparsity"} (proportion zeros) or
#'               \code{"richness"} (number of non-zero features). Default is "sparsity".
#' @param outlier_method Method for detecting outliers: \code{"mad"}, \code{"iqr"}, or
#'                       \code{"both"}. Default is "mad".
#' @param multiplier For MAD method: samples with residuals beyond (median + multiplier * MAD)
#'                   are flagged. Default is 3.
#' @param iqr_multiplier For IQR method: samples with residuals beyond Q1 - multiplier*IQR or
#'                       Q3 + multiplier*IQR are flagged. Default is 1.5 (Tukey's fences).
#' @param direction Which outliers to flag: \code{"high_sparsity"} (more zeros than expected),
#'                  \code{"low_sparsity"} (fewer zeros than expected), or \code{"both"}.
#'                  Default is "high_sparsity" since high sparsity at low depth is concerning.
#' @param verbose Logical. Print summary and recommendations. Default is TRUE.
#'
#' @return A list containing:
#'   \item{sample_metrics}{data.frame with sample_name, depth, sparsity, richness, residual}
#'   \item{outliers}{data.frame of flagged samples with reason}
#'   \item{n_outliers}{Number of outliers detected}
#'   \item{fit_summary}{List with regression coefficients and R-squared}
#'   \item{thresholds}{List with cutoff values used for outlier detection}
#'   \item{recommendation}{Character string with filtering advice}
#'
#' @export
#'
#' @examples
#' # Identify high-sparsity outliers using MAD method
#' # result <- analyze_depth_sparsity(my_table, metric = "sparsity", method = "mad")
#'
#' # Check for outliers in richness vs depth relationship
#' # result <- analyze_depth_sparsity(my_table, metric = "richness", method = "iqr")
#'
#' # Get both types of outliers
#' # result <- analyze_depth_sparsity(my_table, outlier_method = "both", direction = "both")
analyze_depth_sparsity <- function(table, metric = c("sparsity", "richness"),
                                    outlier_method = c("mad", "iqr", "both"),
                                    multiplier = 3, iqr_multiplier = 1.5,
                                    direction = c("high_sparsity", "low_sparsity", "both"),
                                    verbose = TRUE) {
  # Validate inputs
  if (!is.data.frame(table) && !is.matrix(table)) {
    stop("table must be a data.frame or matrix")
  }

  if (ncol(table) < 2) {
    stop("table must have at least one sample column")
  }

  metric <- match.arg(metric)
  outlier_method <- match.arg(outlier_method)
  direction <- match.arg(direction)

  # Extract sample data
  sample_data <- table[, -1, drop = FALSE]
  sample_names <- colnames(sample_data)
  n_features <- nrow(sample_data)

  # Calculate metrics per sample
  depths <- colSums(sample_data)
  richness <- colSums(sample_data > 0)
  n_zeros <- colSums(sample_data == 0)
  sparsity <- n_zeros / n_features

  if (metric == "sparsity") {
    values <- sparsity  # Sparsity as proportion
    value_name <- "sparsity"
  } else {
    values <- richness  # Richness (observed ASVs)
    value_name <- "richness"
  }

  # Create metrics dataframe
  sample_metrics <- data.frame(
    sample_name = sample_names,
    depth = depths,
    richness = richness,
    sparsity = sparsity,
    stringsAsFactors = FALSE
  )

  # Log-transform depth for better linear fit (add small constant to avoid log(0))
  log_depth <- log10(depths + 1)
  log_values <- if (metric == "sparsity") {
    # For sparsity, use logit transformation if values are in (0, 1)
    # Otherwise fall back to raw values
    valid_sparsity <- values > 0 & values < 1
    transformed <- rep(NA, length(values))
    transformed[valid_sparsity] <- log(values[valid_sparsity] / (1 - values[valid_sparsity]))
    # Fill NA with boundary values
    transformed[!valid_sparsity & values == 0] <- -3  # Approx logit(0.05)
    transformed[!valid_sparsity & values == 1] <- 3   # Approx logit(0.95)
    transformed
  } else {
    log10(values + 1)
  }

  # Fit robust linear model (using median-based approach)
  # Sort by depth and find median slope
  ord <- order(log_depth)
  sorted_x <- log_depth[ord]
  sorted_y <- log_values[ord]

  # Use Theil-Sen estimator (median of pairwise slopes) for robustness
  if (length(sorted_x) >= 2) {
    slopes <- outer(sorted_y, sorted_y, "-") / (outer(sorted_x, sorted_x, "-") + .Machine$double.eps)
    # Upper triangle only (avoid self-comparison and duplicates)
    idx <- upper.tri(slopes)
    median_slope <- median(slopes[idx], na.rm = TRUE)

    # Intercept: median of (y - slope * x)
    intercept <- median(sorted_y - median_slope * sorted_x, na.rm = TRUE)
  } else {
    median_slope <- 0
    intercept <- mean(log_values, na.rm = TRUE)
  }

  # Calculate predicted values and residuals
  predicted <- intercept + median_slope * log_depth
  residuals <- log_values - predicted

  sample_metrics$residual <- residuals
  sample_metrics$log_depth <- log_depth
  sample_metrics$value <- values

  # Outlier detection
  outliers_list <- list()

  if (outlier_method %in% c("mad", "both")) {
    med_residual <- median(residuals, na.rm = TRUE)
    mad_residual <- mad(residuals, na.rm = TRUE)

    # Initialize cutoffs
    high_cutoff <- NA
    low_cutoff <- NA

    if (direction == "high_sparsity" || direction == "both") {
      # High sparsity = more zeros than expected = positive residual for sparsity
      # For richness, low richness = negative residual
      if (metric == "sparsity") {
        high_cutoff <- med_residual + multiplier * mad_residual
        high_mask <- residuals > high_cutoff
      } else {
        high_cutoff <- med_residual - multiplier * mad_residual  # Low richness is bad
        high_mask <- residuals < high_cutoff
      }
      high_outliers <- sample_metrics[high_mask, , drop = FALSE]
      if (nrow(high_outliers) > 0) {
        high_outliers$outlier_type <- "high_sparsity_low_richness"
      }
      outliers_list$high <- high_outliers
    }

    if (direction == "low_sparsity" || direction == "both") {
      if (metric == "sparsity") {
        low_cutoff <- med_residual - multiplier * mad_residual
        low_mask <- residuals < low_cutoff
      } else {
        low_cutoff <- med_residual + multiplier * mad_residual  # High richness is good, not an outlier
        low_mask <- residuals > low_cutoff
      }
      low_outliers <- sample_metrics[low_mask, , drop = FALSE]
      if (nrow(low_outliers) > 0) {
        low_outliers$outlier_type <- "low_sparsity_high_richness"
      }
      outliers_list$low <- low_outliers
    }

    attr(sample_metrics, "mad_cutoff_high") <- if (metric == "sparsity") high_cutoff else low_cutoff
    attr(sample_metrics, "mad_cutoff_low") <- if (metric == "sparsity") low_cutoff else high_cutoff
  }

  if (outlier_method %in% c("iqr", "both")) {
    q1 <- quantile(residuals, 0.25, na.rm = TRUE)
    q3 <- quantile(residuals, 0.75, na.rm = TRUE)
    iqr <- q3 - q1

    if (direction == "high_sparsity" || direction == "both") {
      if (metric == "sparsity") {
        iqr_high_cutoff <- q3 + iqr_multiplier * iqr
        iqr_high_mask <- residuals > iqr_high_cutoff
      } else {
        iqr_high_cutoff <- q1 - iqr_multiplier * iqr
        iqr_high_mask <- residuals < iqr_high_cutoff
      }
      iqr_high_outliers <- sample_metrics[iqr_high_mask, , drop = FALSE]
      if (nrow(iqr_high_outliers) > 0) {
        iqr_high_outliers$outlier_type <- "high_sparsity_low_richness_iqr"
      }
      outliers_list$iqr_high <- iqr_high_outliers
    }

    if (direction == "low_sparsity" || direction == "both") {
      if (metric == "sparsity") {
        iqr_low_cutoff <- q1 - iqr_multiplier * iqr
        iqr_low_mask <- residuals < iqr_low_cutoff
      } else {
        iqr_low_cutoff <- q3 + iqr_multiplier * iqr
        iqr_low_mask <- residuals > iqr_low_cutoff
      }
      iqr_low_outliers <- sample_metrics[iqr_low_mask, , drop = FALSE]
      if (nrow(iqr_low_outliers) > 0) {
        iqr_low_outliers$outlier_type <- "low_sparsity_high_richness_iqr"
      }
      outliers_list$iqr_low <- iqr_low_outliers
    }

    attr(sample_metrics, "iqr_q1") <- q1
    attr(sample_metrics, "iqr_q3") <- q3
    attr(sample_metrics, "iqr_range") <- iqr
  }

  # Combine outliers based on method
  if (outlier_method == "mad") {
    all_outliers <- rbind(
      if ("high" %in% names(outliers_list)) outliers_list$high else data.frame(),
      if ("low" %in% names(outliers_list)) outliers_list$low else data.frame()
    )
  } else if (outlier_method == "iqr") {
    all_outliers <- rbind(
      if ("iqr_high" %in% names(outliers_list)) outliers_list$iqr_high else data.frame(),
      if ("iqr_low" %in% names(outliers_list)) outliers_list$iqr_low else data.frame()
    )
  } else {  # both
    all_outliers <- rbind(
      if ("high" %in% names(outliers_list)) outliers_list$high else data.frame(),
      if ("low" %in% names(outliers_list)) outliers_list$low else data.frame(),
      if ("iqr_high" %in% names(outliers_list)) outliers_list$iqr_high else data.frame(),
      if ("iqr_low" %in% names(outliers_list)) outliers_list$iqr_low else data.frame()
    )
  }

  if (nrow(all_outliers) == 0) {
    all_outliers <- data.frame(
      sample_name = character(),
      depth = numeric(),
      richness = numeric(),
      sparsity = numeric(),
      residual = numeric(),
      outlier_type = character(),
      stringsAsFactors = FALSE
    )
  }

  # Calculate fit statistics
  ss_res <- sum(residuals^2, na.rm = TRUE)
  ss_tot <- sum((log_values - mean(log_values, na.rm = TRUE))^2, na.rm = TRUE)
  r_squared <- 1 - (ss_res / (ss_tot + .Machine$double.eps))

  fit_summary <- list(
    slope = median_slope,
    intercept = intercept,
    r_squared = r_squared,
    n_samples = nrow(sample_metrics)
  )

  # Thresholds used
  thresholds <- list(
    method = outlier_method,
    multiplier = multiplier,
    iqr_multiplier = iqr_multiplier,
    direction = direction
  )

  # Generate recommendation
  n_outliers <- nrow(all_outliers)
  pct_outliers <- (n_outliers / nrow(sample_metrics)) * 100

  if (n_outliers > 0) {
    outlier_samples <- paste(all_outliers$sample_name, collapse = ", ")

    if (pct_outliers > 10) {
      recommendation <- sprintf(
        "WARNING: %.1f%% (%d) samples flagged as outliers. Consider reviewing these samples: %s. High sparsity relative to depth may indicate poor sequencing quality.",
        pct_outliers, n_outliers, outlier_samples
      )
    } else {
      recommendation <- sprintf(
        "%d samples (%.1f%%) flagged as outliers: %s. These samples have unusual depth-sparsity relationships.",
        n_outliers, pct_outliers, outlier_samples
      )
    }
  } else {
    recommendation <- sprintf(
      "No outliers detected. All samples follow expected depth-%s relationship.",
      value_name
    )
  }

  if (verbose && n_outliers > 0) {
    message(sprintf("\nDepth-Sparsity Analysis Complete"))
    message(sprintf("  Total samples: %d", nrow(sample_metrics)))
    message(sprintf("  Outliers detected: %d (%.1f%%)", n_outliers, pct_outliers))
    message(sprintf("  Fit R-squared: %.3f", r_squared))
    message(sprintf("  Recommendation: %s", recommendation))
  }

  return(list(
    sample_metrics = sample_metrics,
    outliers = all_outliers,
    n_outliers = n_outliers,
    fit_summary = fit_summary,
    thresholds = thresholds,
    recommendation = recommendation
  ))
}


#' Plot Total Reads vs Observed ASVs with outlier detection
#'
#' Simple visualization showing the relationship between sequencing depth (total reads)
#' and observed richness (number of non-zero ASVs/features). Uses MAD-based outlier
#' detection to flag samples that have unusually low richness for their sequencing depth.
#'
#' This plot helps identify:
#' - Low-depth samples with poor feature discovery (potential undersampling)
#' - High-depth samples with unexpectedly few features (potential technical artifacts)
#' - Samples that deviate from the expected depth-richness curve
#'
#' @param table A feature table (data.frame or matrix) with features as rows and samples as columns.
#'              First column should be feature IDs, remaining columns are sample counts.
#' @param mad_multiplier Multiplier for MAD-based outlier detection. Samples with residuals
#'                       beyond (median + multiplier * MAD) are flagged. Default is 3.
#' @param main Main title for the plot. Default is "Total Reads vs Observed ASVs".
#' @param color Theme color. Default is "blue". Options: "blue", "green", "red", "purple"
#' @param show_labels Logical. Show sample names for outliers. Default is TRUE if <= 10 outliers.
#'
#' @return A list containing:
#'   \item{plot}{ggplot object}
#'   \item{outliers}{data.frame of flagged outlier samples}
#'   \item{metrics}{data.frame with all sample metrics}
#'
#' @export
#'
#' @examples
#' # Simple plot with automatic outlier detection
#' # result <- plot_reads_vs_asvs(my_table)
#' # print(result$plot)
#'
#' # Adjust sensitivity
#' # result <- plot_reads_vs_asvs(my_table, mad_multiplier = 2.5)
#'
#' # Access outlier information
#' # outlier_samples <- result$outliers$sample_name
plot_reads_vs_asvs <- function(table, mad_multiplier = 3,
                                main = "Total Reads vs Observed ASVs",
                                color = "blue", show_labels = NULL) {
  # Check for ggplot2
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting. Please install it.")
  }

  library(ggplot2)

  # Validate input
  if (!is.data.frame(table) && !is.matrix(table)) {
    stop("table must be a data.frame or matrix")
  }

  if (ncol(table) < 2) {
    stop("table must have at least one sample column")
  }

  # Calculate per-sample metrics
  sample_data <- table[, -1, drop = FALSE]
  sample_names <- colnames(sample_data)

  depths <- colSums(sample_data)
  observed_asvs <- colSums(sample_data > 0)

  metrics <- data.frame(
    sample_name = sample_names,
    total_reads = depths,
    observed_asvs = observed_asvs,
    stringsAsFactors = FALSE
  )

  # Log-transform for linear relationship
  log_reads <- log10(depths + 1)
  log_asvs <- log10(observed_asvs + 1)

  # Theil-Sen robust slope estimator
  ord <- order(log_reads)
  sorted_x <- log_reads[ord]
  sorted_y <- log_asvs[ord]

  if (length(sorted_x) >= 2) {
    slopes <- outer(sorted_y, sorted_y, "-") / (outer(sorted_x, sorted_x, "-") + .Machine$double.eps)
    idx <- upper.tri(slopes)
    median_slope <- median(slopes[idx], na.rm = TRUE)
    intercept <- median(sorted_y - median_slope * sorted_x, na.rm = TRUE)
  } else {
    median_slope <- 1
    intercept <- 0
  }

  # Calculate residuals
  predicted <- intercept + median_slope * log_reads
  residuals <- log_asvs - predicted

  metrics$residual <- residuals
  metrics$log_reads <- log_reads
  metrics$log_asvs <- log_asvs

  # MAD-based outlier detection
  med_residual <- median(residuals, na.rm = TRUE)
  mad_residual <- mad(residuals, na.rm = TRUE)

  # Flag samples with unusually LOW richness (negative residuals beyond cutoff)
  low_richness_cutoff <- med_residual - mad_multiplier * mad_residual
  outlier_mask <- residuals < low_richness_cutoff

  metrics$outlier <- outlier_mask
  metrics$outlier_type <- ifelse(outlier_mask, "low_richness", "normal")

  outliers <- metrics[outlier_mask, , drop = FALSE]

  # Auto-determine whether to show labels
  if (is.null(show_labels)) {
    show_labels <- nrow(outliers) <= 10
  }

  # Color mapping
  color_map <- list(
    blue = "#2E86AB",
    green = "#239B56",
    red = "#C0392B",
    purple = "#884EA0"
  )
  primary_color <- color_map[[color]]
  if (is.null(primary_color)) primary_color <- color_map[["blue"]]

  base_theme <- theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      panel.grid.minor = element_line(color = "grey90")
    )

  # Create plot
  p <- ggplot(metrics, aes(x = total_reads, y = observed_asvs)) +
    geom_point(aes(color = outlier_type), alpha = 0.6, size = 2) +
    scale_color_manual(
      name = "",
      values = c("normal" = "grey50", "low_richness" = "#C0392B"),
      limits = c("normal", "low_richness")
    ) +
    geom_abline(
      slope = 10^(intercept) * median_slope,
      intercept = 0,
      color = primary_color,
      linewidth = 1,
      linetype = "dashed"
    ) +
    labs(
      title = main,
      subtitle = sprintf(
        "MAD multiplier = %.1f | %d outliers flagged (low richness for depth)",
        mad_multiplier, sum(outlier_mask)
      ),
      x = "Total Reads (sequencing depth)",
      y = "Observed ASVs (richness)"
    ) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    base_theme

  # Add sample name labels for outliers
  if (show_labels && nrow(outliers) > 0) {
    p <- p +
      geom_text(
        data = outliers,
        aes(x = total_reads * 1.03, y = observed_asvs * 1.02, label = sample_name),
        size = 3, hjust = 0, vjust = 0, angle = 0, color = "#C0392B"
      )
  }

  return(list(
    plot = p,
    outliers = outliers,
    metrics = metrics,
    cutoff = list(
      method = "MAD",
      multiplier = mad_multiplier,
      residual_threshold = low_richness_cutoff
    )
  ))
}


#' Plot depth-sparsity relationship with outlier detection
#'
#' Creates visualization of the relationship between sequencing depth and
#' sparsity/richness, with outliers highlighted.
#'
#' @param analysis_result Output from \code{\link{analyze_depth_sparsity}}
#' @param main Main title for the plot
#' @param color Theme color. Default is "blue". Options: "blue", "green", "red", "purple"
#'
#' @return A ggplot object
#'
#' @export
#'
#' @examples
#' # result <- analyze_depth_sparsity(my_table)
#' # plot_depth_sparsity(result)
plot_depth_sparsity <- function(analysis_result, main = "Depth-Sparsity Relationship",
                                 color = "blue") {
  # Check for ggplot2
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting. Please install it.")
  }

  library(ggplot2)

  metrics_df <- analysis_result$sample_metrics
  outliers_df <- analysis_result$outliers
  fit <- analysis_result$fit_summary

  # Color mapping
  color_map <- list(
    blue = "#2E86AB",
    green = "#239B56",
    red = "#C0392B",
    purple = "#884EA0"
  )
  primary_color <- color_map[[color]]
  if (is.null(primary_color)) primary_color <- color_map[["blue"]]

  base_theme <- theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      panel.grid.minor = element_line(color = "grey90")
    )

  # Determine y-axis label based on outlier type
  has_high_sparsity <- any(grepl("high_sparsity", outliers_df$outlier_type, ignore.case = TRUE))

  p <- ggplot(metrics_df, aes(x = depth, y = sparsity * 100)) +
    geom_point(aes(color = "All Samples"), alpha = 0.6, size = 2) +
    geom_line(aes(y = exp(fit$intercept + fit$slope * log10(depth + 1)) /
                        (1 + exp(fit$intercept + fit$slope * log10(depth + 1))) * 100,
                  color = "Fitted Line"), linewidth = 1) +
    geom_point(data = outliers_df,
               aes(x = depth, y = sparsity * 100, color = "Outliers"),
               size = 3, shape = 17) +
    scale_color_manual(
      name = "",
      values = c(
        "All Samples" = "grey50",
        "Fitted Line" = primary_color,
        "Outliers" = "#C0392B"
      )
    ) +
    labs(
      title = main,
      subtitle = sprintf(
        "R² = %.3f | Slope = %.3f | %d outliers detected",
        fit$r_squared, fit$slope, analysis_result$n_outliers
      ),
      x = "Sequencing Depth (reads)",
      y = "Sparsity (% zeros)"
    ) +
    scale_x_log10() +
    base_theme

  # Add text labels for extreme outliers
  if (nrow(outliers_df) > 0 && nrow(outliers_df) <= 10) {
    p <- p +
      geom_text(
        data = outliers_df,
        aes(x = depth * 1.05, y = sparsity * 100 * 1.02, label = sample_name),
        size = 3, hjust = 0, vjust = 0, angle = 0
      )
  }

  return(p)
}


#' Convenience function: filter samples based on depth-sparsity outliers
#'
#' Combines \code{\link{analyze_depth_sparsity}} with filtering to remove
#' flagged outlier samples in one step.
#'
#' @param table A feature table (data.frame or matrix)
#' @param metric What to analyze: \code{"sparsity"} or \code{"richness"}
#' @param outlier_method Method: \code{"mad"}, \code{"iqr"}, or \code{"both"}
#' @param multiplier MAD multiplier for outlier detection
#' @param keep_outliers Logical. If TRUE, keep outliers; if FALSE (default), remove them.
#'
#' @return Filtered feature table (same structure as input)
#'
#' @export
#'
#' @examples
#' # Remove high-sparsity outlier samples
#' # cleaned_table <- filter_depth_sparsity_outliers(my_table)
#'
#' # Keep only samples that pass the depth-sparsity check
#' # cleaned_table <- filter_depth_sparsity_outliers(my_table, keep_outliers = FALSE)
filter_depth_sparsity_outliers <- function(table, metric = c("sparsity", "richness"),
                                            outlier_method = c("mad", "iqr", "both"),
                                            multiplier = 3, keep_outliers = FALSE) {
  metric <- match.arg(metric)
  outlier_method <- match.arg(outlier_method)

  # Run analysis
  result <- analyze_depth_sparsity(
    table, metric = metric, outlier_method = outlier_method,
    multiplier = multiplier, verbose = FALSE
  )

  # Get outlier sample names
  outlier_names <- result$outliers$sample_name

  # Determine which samples to keep
  if (keep_outliers) {
    keep_samples <- outlier_names
  } else {
    all_samples <- colnames(table[, -1, drop = FALSE])
    keep_samples <- setdiff(all_samples, outlier_names)
  }

  # Filter table
  if (length(keep_samples) == 0) {
    warning("No samples would remain after filtering. Returning original table.")
    return(table)
  }

  result_table <- table[, c(TRUE, keep_samples), drop = FALSE]

  # Preserve column names
  colnames(result_table) <- colnames(table)[c(TRUE, keep_samples)]

  # Attach analysis results as attributes
  attr(result_table, "outlier_analysis") <- result
  attr(result_table, "n_removed") <- if (keep_outliers) length(setdiff(colnames(table)[-1], keep_samples)) else length(outlier_names)
  attr(result_table, "n_retained") <- ncol(result_table) - 1

  return(result_table)
}
