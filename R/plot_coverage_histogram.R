#' Plot histogram of sample coverage using ggplot2
#'
#' Creates a histogram of total reads per sample with dynamic bin sizing.
#' Optionally highlights samples above/below a threshold.
#'
#' @param table A feature table (data.frame or matrix) with features as rows and samples as columns.
#'              First column should be feature IDs, remaining columns are sample counts.
#' @param threshold Optional coverage threshold. Bars below this value are colored red (to remove),
#'                  bars at or above are colored green (to keep).
#' @param main Title for the plot. Default is "Sample Coverage Distribution".
#' @param xlab Label for x-axis. Default is "Total Reads per Sample".
#' @param col_above Color for bars at or above threshold (default: "steelblue").
#' @param col_below Color for bars below threshold (default: "tomato").
#' @param plot_dir Optional directory path to save the plot. If NULL, plot is displayed only.
#' @param prefix Prefix for saved plot filename. Default is "coverage_histogram".
#' @param width Plot width in inches for ggsave. Default 8.
#' @param height Plot height in inches for ggsave. Default 6.
#' @param dpi DPI for saved raster output. Default 300.
#'
#' @return Returns a ggplot object invisibly. If plot_dir is provided, also saves the plot.
#'
#' @export
#'
#' @examples
#' data(example_feature_table)
#' plot_coverage_histogram(example_feature_table)
plot_coverage_histogram <- function(table, threshold = NULL,
                                     main = "Sample Coverage Distribution",
                                     xlab = "Total Reads per Sample",
                                     col_above = "steelblue",
                                     col_below = "tomato",
                                     plot_dir = NULL,
                                     prefix = "coverage_histogram",
                                     width = 8, height = 6, dpi = 300) {
  # Calculate total reads per sample
  coverage_vec <- colSums(table[, -1, drop = FALSE])

  # Create data frame for ggplot
  df <- data.frame(coverage = coverage_vec)

  # Dynamically determine number of bins using Sturges' rule (default)
  n_bins <- ceiling(log2(length(coverage_vec)) + 1)

  # Create histogram object without plotting to get breaks
  hist_obj <- hist(coverage_vec, breaks = n_bins, plot = FALSE)

  # Adjust bin width if all values are the same
  if (diff(range(hist_obj$breaks)) == 0) {
    hist_obj$breaks <- c(min(coverage_vec) - 0.5, max(coverage_vec) + 0.5)
    hist_obj$counts <- length(coverage_vec)
    hist_obj$mids <- mean(coverage_vec)
  }

  # Create data frame with bin information
  bin_data <- data.frame(
    mid = hist_obj$mids,
    count = hist_obj$counts,
    break_min = head(hist_obj$breaks, -1),
    break_max = tail(hist_obj$breaks, -1)
  )

  # Add color category based on threshold
  if (!is.null(threshold)) {
    bin_data$category <- ifelse(bin_data$mid >= threshold, "Keep", "Remove")
    bin_data$color_label <- ifelse(bin_data$mid >= threshold,
                                    paste0("Keep (>= ", threshold, ")"),
                                    paste0("Remove (< ", threshold, ")"))
    plot_colors <- c("Keep" = col_above, "Remove" = col_below)
  } else {
    bin_data$category <- "Coverage"
    bin_data$color_label <- "Coverage"
    plot_colors <- c("Coverage" = col_above)
  }

  # Build ggplot
  p <- ggplot2::ggplot(bin_data, ggplot2::aes(x = mid, y = count, fill = category)) +
    ggplot2::geom_col(width = diff(hist_obj$breaks)[1],
                       color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = plot_colors, na.value = col_above) +
    ggplot2::labs(title = main, x = xlab, y = "Number of Samples", fill = NULL) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  # Add legend positioning
  if (!is.null(threshold)) {
    p <- p + ggplot2::theme(legend.position = "top")
  }

  # Save plot if directory specified
  if (!is.null(plot_dir)) {
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE)
    }
    ggplot2::ggsave(
      filename = file.path(plot_dir, paste0(prefix, ".png")),
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
  }

  # Return invisible list with coverage data and breaks (backward compatible)
  invisible(list(
    coverage = coverage_vec,
    breaks = hist_obj$breaks,
    counts = hist_obj$counts,
    plot = p
  ))
}
