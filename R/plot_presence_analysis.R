#' Plot presence frequency analysis for features and samples using ggplot2
#'
#' Creates histograms showing:
#' - In how many samples each feature is present (feature prevalence)
#' - How many features are present in each sample (sample richness)
#'
#' @param table A feature table (data.frame or matrix) with features as rows and samples as columns.
#'              First column should be feature IDs, remaining columns are sample counts.
#' @param threshold Minimum abundance for a feature to be considered "present" in a sample.
#'                  Default is 1 (any non-zero count).
#' @param save_dir Optional directory to save plots. If NULL, plots are displayed only.
#' @param prefix Prefix for saved plot filenames. Default is "presence_analysis".
#' @param main_main Main title for the combined plot. Default is "Presence Frequency Analysis".
#' @param width Plot width in inches for ggsave. Default 10.
#' @param height Plot height in inches for ggsave. Default 6.
#' @param dpi DPI for saved raster output. Default 300.
#' @param prev_bins Number of bins for feature prevalence histogram. Default is 20.
#' @param rich_bins Number of bins for sample richness histogram. Default is 20.
#' @param table_filtered Optional filtered feature table for comparison. If provided,
#'                       shows original vs filtered side by side.
#' @param prefix_filtered Prefix for filtered plot files. Default is "filtered".
#'
#' @return Returns a list containing:
#' \describe{
#'   \item{feature_prevalence}{Named vector: number of samples each feature appears in (original)}
#'   \item{sample_richness}{Named vector: number of features in each sample (original)}
#'   \item{feature_prevalence_filtered}{Named vector: filtered feature prevalence (if table_filtered provided)}
#'   \item{sample_richness_filtered}{Named vector: filtered sample richness (if table_filtered provided)}
#'   \item{mean_feature_prevalence}{Mean number of samples per feature (original)}
#'   \item{mean_sample_richness}{Mean number of features per sample (original)}
#'   \item{plots}{List containing the ggplot objects}
#' }
#'
#' @export
#'
#' @examples
#' data(example_feature_table)
#' plot_presence_analysis(example_feature_table)
plot_presence_analysis <- function(table, threshold = 1,
                                    save_dir = NULL, prefix = "presence_analysis",
                                    main_main = "Presence Frequency Analysis",
                                    width = 10, height = 6, dpi = 300,
                                    prev_bins = 20, rich_bins = 20,
                                    table_filtered = NULL, prefix_filtered = "filtered") {
  # Extract abundance data
  feature_ids <- table[, 1]
  abundances <- as.matrix(table[, -1, drop = FALSE])

  # Create presence/absence matrix (1 = present, 0 = absent)
  presence <- (abundances >= threshold) * 1

  # Calculate feature prevalence (in how many samples each feature appears)
  # colSums counts non-zero values per column (sample) -> but we want per feature across samples
  # Actually: presence matrix has features as rows, samples as columns
  # So colSums gives: for each sample, how many features are present = sample richness
  # And rowSums gives: for each feature, how many samples it appears in = feature prevalence
  feature_prevalence <- rowSums(presence)  # Each feature: count samples where it's present
  names(feature_prevalence) <- feature_ids

  # Calculate sample richness (how many features in each sample)
  sample_richness <- colSums(presence)  # Each sample: count features that are present
  names(sample_richness) <- colnames(abundances)

  # Calculate mean and median for annotation
  mean_prev <- mean(feature_prevalence)
  median_prev <- median(feature_prevalence)
  mean_rich <- mean(sample_richness)
  median_rich <- median(sample_richness)

  # Create data frames for histogram plotting
  prev_df <- data.frame(
    value = feature_prevalence
  )

  # Sample richness frequency
  rich_df <- data.frame(
    value = sample_richness
  )

  # Create output directory if specified
  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
  }

  # Helper function to create prevalence plot with comparison option
  create_prevalence_plot <- function(prev_df, mean_val, median_val, bins, title_suffix = "") {
    ggplot2::ggplot(prev_df, ggplot2::aes(x = value, fill = ..fill.color..)) +
      ggplot2::geom_histogram(bins = bins, color = "white", linewidth = 0.5) +
      ggplot2::geom_vline(
        data = data.frame(x = c(mean_val, median_val), type = c("Mean", "Median")),
        ggplot2::aes(xintercept = x, color = type, linetype = type),
        linewidth = 1, show.legend = TRUE
      ) +
      ggplot2::labs(
        title = paste0("Feature Prevalence", title_suffix, "\n(In How Many Samples)"),
        x = "Number of Samples Feature Appears In",
        y = "Number of Features",
        color = "Statistics",
        linetype = "Statistics"
      ) +
      ggplot2::scale_color_manual(
        values = c("Mean" = "black", "Median" = "gray40"),
        guide = ggplot2::guide_legend(order = 1)
      ) +
      ggplot2::scale_linetype_manual(
        values = c("Mean" = "dashed", "Median" = "dotdash"),
        guide = ggplot2::guide_legend(order = 2)
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "top",
        legend.box = "vertical"
      )
  }

  # Helper function to create richness plot with comparison option
  create_richness_plot <- function(rich_df, mean_val, median_val, bins, title_suffix = "") {
    ggplot2::ggplot(rich_df, ggplot2::aes(x = value, fill = ..fill.color..)) +
      ggplot2::geom_histogram(bins = bins, color = "white", linewidth = 0.5) +
      ggplot2::geom_vline(
        data = data.frame(x = c(mean_val, median_val), type = c("Mean", "Median")),
        ggplot2::aes(xintercept = x, color = type, linetype = type),
        linewidth = 1, show.legend = TRUE
      ) +
      ggplot2::labs(
        title = paste0("Sample Richness", title_suffix, "\n(How Many Features Per Sample)"),
        x = "Number of Features in Sample",
        y = "Number of Samples",
        color = "Statistics",
        linetype = "Statistics"
      ) +
      ggplot2::scale_color_manual(
        values = c("Mean" = "black", "Median" = "gray40"),
        guide = ggplot2::guide_legend(order = 1)
      ) +
      ggplot2::scale_linetype_manual(
        values = c("Mean" = "dashed", "Median" = "dotdash"),
        guide = ggplot2::guide_legend(order = 2)
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "top",
        legend.box = "vertical"
      )
  }

  plots_list <- list()

  if (!is.null(table_filtered)) {
    # Process filtered table
    feature_ids_filt <- table_filtered[, 1]
    abundances_filt <- as.matrix(table_filtered[, -1, drop = FALSE])
    presence_filt <- (abundances_filt >= threshold) * 1

    feature_prevalence_filt <- rowSums(presence_filt)
    names(feature_prevalence_filt) <- feature_ids_filt
    sample_richness_filt <- colSums(presence_filt)
    names(sample_richness_filt) <- colnames(abundances_filt)

    mean_prev_filt <- mean(feature_prevalence_filt)
    median_prev_filt <- median(feature_prevalence_filt)
    mean_rich_filt <- mean(sample_richness_filt)
    median_rich_filt <- median(sample_richness_filt)

    prev_df_filt <- data.frame(value = feature_prevalence_filt, group = "Filtered")
    rich_df_filt <- data.frame(value = sample_richness_filt, group = "Filtered")

    # Combined data for overlay plots
    prev_df_combined <- data.frame(
      value = c(feature_prevalence, feature_prevalence_filt),
      group = rep(c("Original", "Filtered"), c(length(feature_prevalence), length(feature_prevalence_filt)))
    )
    rich_df_combined <- data.frame(
      value = c(sample_richness, sample_richness_filt),
      group = rep(c("Original", "Filtered"), c(length(sample_richness), length(sample_richness_filt)))
    )

    # Calculate common x-axis limits for both plots
    max_prev_x <- max(max(feature_prevalence), max(feature_prevalence_filt))
    max_rich_x <- max(max(sample_richness), max(sample_richness_filt))

    # Create data frames for statistics with proper facet mapping
    prev_stats_df <- data.frame(
      x = c(mean_prev, median_prev, mean_prev_filt, median_prev_filt),
      stat = rep(c("Mean", "Median"), 2),
      group = rep(c("Original", "Filtered"), each = 2)
    )

    rich_stats_df <- data.frame(
      x = c(mean_rich, median_rich, mean_rich_filt, median_rich_filt),
      stat = rep(c("Mean", "Median"), 2),
      group = rep(c("Original", "Filtered"), each = 2)
    )

    # Plot 1: Feature Prevalence comparison (side by side using facets)
    p1 <- ggplot2::ggplot(prev_df_combined, ggplot2::aes(x = value, fill = group)) +
      ggplot2::geom_histogram(bins = prev_bins, color = "white", linewidth = 0.5) +
      ggplot2::facet_wrap(~group, ncol = 2) +
      ggplot2::geom_vline(
        data = prev_stats_df,
        ggplot2::aes(xintercept = x, linetype = stat),
        color = "black",
        linewidth = 1
      ) +
      ggplot2::labs(
        title = "Feature Prevalence Comparison",
        x = "Number of Samples Feature Appears In",
        y = "Number of Features",
        fill = "Dataset",
        linetype = "Statistic"
      ) +
      ggplot2::scale_linetype_manual(values = c("Mean" = "dashed", "Median" = "dotdash")) +
      ggplot2::scale_fill_manual(values = c("Original" = "steelblue", "Filtered" = "coral")) +
      ggplot2::coord_cartesian(xlim = c(0, max_prev_x)) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "top",
        strip.background = ggplot2::element_rect(fill = "gray90"),
        strip.text = ggplot2::element_text(face = "bold")
      )

    # Plot 2: Sample Richness comparison (side by side using facets)
    p2 <- ggplot2::ggplot(rich_df_combined, ggplot2::aes(x = value, fill = group)) +
      ggplot2::geom_histogram(bins = rich_bins, color = "white", linewidth = 0.5) +
      ggplot2::facet_wrap(~group, ncol = 2) +
      ggplot2::geom_vline(
        data = rich_stats_df,
        ggplot2::aes(xintercept = x, linetype = stat),
        color = "black",
        linewidth = 1
      ) +
      ggplot2::labs(
        title = "Sample Richness Comparison",
        x = "Number of Features in Sample",
        y = "Number of Samples",
        fill = "Dataset",
        linetype = "Statistic"
      ) +
      ggplot2::scale_linetype_manual(values = c("Mean" = "dashed", "Median" = "dotdash")) +
      ggplot2::scale_fill_manual(values = c("Original" = "steelblue", "Filtered" = "coral")) +
      ggplot2::coord_cartesian(xlim = c(0, max_rich_x)) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "top",
        strip.background = ggplot2::element_rect(fill = "gray90"),
        strip.text = ggplot2::element_text(face = "bold")
      )

    plots_list <- list(prevalence = p1, richness = p2)

    # Save combined plot if directory specified
    if (!is.null(save_dir)) {
      if (requireNamespace("patchwork", quietly = TRUE)) {
        combined_plot <- patchwork::wrap_plots(p1, p2, ncol = 2, widths = c(1, 1))
        combined_plot <- combined_plot + patchwork::plot_annotation(
          title = main_main,
          tag_levels = "A",
          tag_prefix = "",
          theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold"))
        )
      } else {
        # Fallback: save only the first plot
        combined_plot <- p1
        warning("patchwork not available. Saving only prevalence plot. Install patchwork for combined plots.")
      }

      ggplot2::ggsave(
        filename = file.path(save_dir, paste0(prefix, "_presence_frequency_comparison.png")),
        plot = combined_plot,
        width = width * 1.2,
        height = height,
        dpi = dpi
      )
    } else {
      if (requireNamespace("patchwork", quietly = TRUE)) {
        combined_plot <- patchwork::wrap_plots(p1, p2, ncol = 2)
        print(combined_plot)
      } else {
        print(p1)
        print(p2)
      }
    }

    # Return stats for both original and filtered
    return(list(
      feature_prevalence = feature_prevalence,
      sample_richness = sample_richness,
      feature_prevalence_filtered = feature_prevalence_filt,
      sample_richness_filtered = sample_richness_filt,
      mean_feature_prevalence = mean_prev,
      mean_sample_richness = mean_rich,
      mean_feature_prevalence_filtered = mean_prev_filt,
      mean_sample_richness_filtered = mean_rich_filt,
      min_feature_prevalence = min(feature_prevalence),
      max_feature_prevalence = max(feature_prevalence),
      min_sample_richness = min(sample_richness),
      max_sample_richness = max(sample_richness),
      min_feature_prevalence_filtered = min(feature_prevalence_filt),
      max_feature_prevalence_filtered = max(feature_prevalence_filt),
      min_sample_richness_filtered = min(sample_richness_filt),
      max_sample_richness_filtered = max(sample_richness_filt),
      n_features = length(feature_ids),
      n_samples = ncol(abundances),
      n_features_filtered = length(feature_ids_filt),
      n_samples_filtered = ncol(abundances_filt),
      plots = plots_list
    ))
  } else {
    # Original single-table behavior
    p1 <- ggplot2::ggplot(prev_df, ggplot2::aes(x = value)) +
      ggplot2::geom_histogram(bins = prev_bins, fill = "steelblue", color = "white", linewidth = 0.5) +
      ggplot2::geom_vline(
        data = data.frame(x = c(mean_prev, median_prev), type = c("Mean", "Median")),
        ggplot2::aes(xintercept = x, color = type, linetype = type),
        linewidth = 1
      ) +
      ggplot2::labs(
        title = "Feature Prevalence\n(In How Many Samples)",
        x = "Number of Samples Feature Appears In",
        y = "Number of Features"
      ) +
      ggplot2::scale_color_manual(
        name = "Statistics",
        values = c("Mean" = "red", "Median" = "orange")
      ) +
      ggplot2::scale_linetype_manual(
        name = "Statistics",
        values = c("Mean" = "dashed", "Median" = "dotdash")
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "top"
      )

    p2 <- ggplot2::ggplot(rich_df, ggplot2::aes(x = value)) +
      ggplot2::geom_histogram(bins = rich_bins, fill = "darkgreen", color = "white", linewidth = 0.5) +
      ggplot2::geom_vline(
        data = data.frame(x = c(mean_rich, median_rich), type = c("Mean", "Median")),
        ggplot2::aes(xintercept = x, color = type, linetype = type),
        linewidth = 1
      ) +
      ggplot2::labs(
        title = "Sample Richness\n(How Many Features Per Sample)",
        x = "Number of Features in Sample",
        y = "Number of Samples"
      ) +
      ggplot2::scale_color_manual(
        name = "Statistics",
        values = c("Mean" = "red", "Median" = "orange")
      ) +
      ggplot2::scale_linetype_manual(
        name = "Statistics",
        values = c("Mean" = "dashed", "Median" = "dotdash")
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "top"
      )

    plots_list <- list(prevalence = p1, richness = p2)

    # Save combined plot if directory specified
    if (!is.null(save_dir)) {
      if (requireNamespace("patchwork", quietly = TRUE)) {
        combined_plot <- patchwork::wrap_plots(p1, p2, ncol = 2, widths = c(1, 1))
        combined_plot <- combined_plot + patchwork::plot_annotation(
          title = main_main,
          tag_levels = "A",
          tag_prefix = "",
          theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold"))
        )
      } else {
        # Fallback: save only the first plot
        combined_plot <- p1
        warning("patchwork not available. Saving only prevalence plot. Install patchwork for combined plots.")
      }

      ggplot2::ggsave(
        filename = file.path(save_dir, paste0(prefix, "_presence_frequency.png")),
        plot = combined_plot,
        width = width,
        height = height,
        dpi = dpi
      )
    } else {
      if (requireNamespace("patchwork", quietly = TRUE)) {
        combined_plot <- patchwork::wrap_plots(p1, p2, ncol = 2)
        print(combined_plot)
      } else {
        print(p1)
        print(p2)
      }
    }

    return(list(
      feature_prevalence = feature_prevalence,
      sample_richness = sample_richness,
      mean_feature_prevalence = mean_prev,
      mean_sample_richness = mean_rich,
      min_feature_prevalence = min(feature_prevalence),
      max_feature_prevalence = max(feature_prevalence),
      min_sample_richness = min(sample_richness),
      max_sample_richness = max(sample_richness),
      n_features = length(feature_ids),
      n_samples = ncol(abundances),
      plots = plots_list
    ))
  }
}
