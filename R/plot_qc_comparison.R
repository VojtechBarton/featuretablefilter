#' Create visualization comparison of original vs filtered feature tables using ggplot2
#'
#' Generates a comprehensive set of plots to compare original and filtered tables.
#' Includes coverage distributions, sparsity histograms, retention barplots, top features barplots, and heatmaps.
#'
#' Plots generated:
#' \itemize{
#'   \item Sample coverage distribution (log10 scale, faceted)
#'   \item Feature abundance distribution (log10 scale, faceted)
#'   \item Sample sparsity histogram (per-sample sparsity distribution, faceted with mean lines)
#'   \item Retention rates barplot (read, feature, sample retention)
#'   \item Top N features stacked barplot (relative abundance comparison)
#'   \item Heatmap of top variable features (original and filtered separately)
#' }
#'
#' @param original_table Original feature table (data.frame) before filtering.
#' @param filtered_table Filtered feature table (data.frame) after filtering.
#' @param plot_dir Optional directory path to save plots. If NULL, plots are displayed only.
#' @param prefix Prefix for saved plot filenames. Default is "qc_comparison".
#' @param top_n Number of top features to show in barplots. Default is 10.
#' @param heatmap_top Number of top variable features to show in heatmap. Default is 50.
#' @param width Plot width in inches for ggsave. Default 10.
#' @param height Plot height in inches for ggsave. Default 8.
#' @param dpi DPI for saved raster output. Default 300.
#'
#' @return A list containing:
#'   \item{plots}{Named list of ggplot objects}
#'   \item{paths}{Named list of file paths (if plot_dir provided)}
#'
#' @export
#'
#' @examples
#' # Display plots
#' # plot_qc_comparison(original_table, filtered_table)
#'
#' # Save plots to directory
#' # plot_qc_comparison(original_table, filtered_table, plot_dir = "qc_plots", prefix = "my_analysis")
plot_qc_comparison <- function(original_table, filtered_table,
                                plot_dir = NULL, prefix = "qc_comparison",
                                top_n = 10, heatmap_top = 50,
                                width = 10, height = 8, dpi = 300) {
  # Create output directory if specified
  if (!is.null(plot_dir)) {
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE)
    }
  }

  # Extract data
  orig_abund <- as.matrix(original_table[, -1, drop = FALSE])
  filt_abund <- as.matrix(filtered_table[, -1, drop = FALSE])
  orig_features <- original_table[, 1]
  filt_features <- filtered_table[, 1]

  # Calculate coverages
  orig_coverage <- colSums(orig_abund)
  filt_coverage <- colSums(filt_abund)

  # Calculate feature abundances (total reads per feature across all samples)
  orig_feature_abund <- rowSums(orig_abund)
  filt_feature_abund <- rowSums(filt_abund)

  # === Plot 1: Sample Coverage Distribution ===
  cov_data <- data.frame(
    coverage = c(log10(orig_coverage + 1), log10(filt_coverage + 1)),
    table_type = rep(c("Original", "Filtered"), times = c(length(orig_coverage), length(filt_coverage)))
  )

  p1 <- ggplot2::ggplot(cov_data, ggplot2::aes(x = coverage, fill = table_type)) +
    ggplot2::geom_histogram(bins = 30, color = "white", linewidth = 0.3, alpha = 0.8) +
    ggplot2::facet_wrap(~table_type, ncol = 2) +
    ggplot2::labs(
      title = paste0(prefix, "_sample_coverage_distribution"),
      x = bquote(log[10]~(Coverage~+~1)),
      y = "Number of Samples",
      fill = NULL
    ) +
    ggplot2::scale_fill_manual(values = c("Original" = "steelblue", "Filtered" = "darkgreen")) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(plot_dir)) {
    ggplot2::ggsave(
      filename = file.path(plot_dir, paste0(prefix, "_sample_coverage_distribution.png")),
      plot = p1,
      width = width * 1.2,
      height = height,
      dpi = dpi
    )
  }

  # === Plot 2: Feature Abundance Distribution ===
  feat_data <- data.frame(
    abundance = c(log10(orig_feature_abund + 1), log10(filt_feature_abund + 1)),
    table_type = rep(c("Original", "Filtered"), times = c(length(orig_feature_abund), length(filt_feature_abund)))
  )

  p_feat <- ggplot2::ggplot(feat_data, ggplot2::aes(x = abundance, fill = table_type)) +
    ggplot2::geom_histogram(bins = 30, color = "white", linewidth = 0.3, alpha = 0.8) +
    ggplot2::facet_wrap(~table_type, ncol = 2) +
    ggplot2::labs(
      title = paste0(prefix, "_feature_abundance_distribution"),
      x = bquote(log[10]~(Abundance~+~1)),
      y = "Number of Features",
      fill = NULL
    ) +
    ggplot2::scale_fill_manual(values = c("Original" = "steelblue", "Filtered" = "darkgreen")) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(plot_dir)) {
    ggplot2::ggsave(
      filename = file.path(plot_dir, paste0(prefix, "_feature_abundance_distribution.png")),
      plot = p_feat,
      width = width * 1.2,
      height = height,
      dpi = dpi
    )
  }

  # === Plot 2: Sparsity Comparison (per sample histogram) ===
  # Calculate sparsity per sample (proportion of zeros in each sample)
  orig_sparsity_per_sample <- colMeans(orig_abund == 0) * 100
  filt_sparsity_per_sample <- colMeans(filt_abund == 0) * 100

  sparsity_data <- data.frame(
    sparsity = c(orig_sparsity_per_sample, filt_sparsity_per_sample),
    table_type = rep(c("Original", "Filtered"), times = c(length(orig_sparsity_per_sample), length(filt_sparsity_per_sample)))
  )

  # Calculate overall means for vertical lines
  orig_mean_sparsity <- mean(orig_sparsity_per_sample)
  filt_mean_sparsity <- mean(filt_sparsity_per_sample)

  sparsity_stats <- data.frame(
    table_type = c("Original", "Filtered"),
    mean_sparsity = c(orig_mean_sparsity, filt_mean_sparsity)
  )

  p2 <- ggplot2::ggplot(sparsity_data, ggplot2::aes(x = sparsity, fill = table_type)) +
    ggplot2::geom_histogram(bins = 20, color = "white", linewidth = 0.5, alpha = 0.8) +
    ggplot2::facet_wrap(~table_type, ncol = 2) +
    ggplot2::geom_vline(
      data = sparsity_stats,
      ggplot2::aes(xintercept = mean_sparsity, color = table_type),
      linewidth = 1, linetype = "dashed"
    ) +
    ggplot2::labs(
      title = "Sample Sparsity Distribution\n(Percentage of Zero Values Per Sample)",
      x = "Sparsity (%)",
      y = "Number of Samples",
      fill = "Dataset",
      color = "Mean"
    ) +
    ggplot2::scale_fill_manual(values = c("Original" = "steelblue", "Filtered" = "darkgreen")) +
    ggplot2::scale_color_manual(values = c("Original" = "darkblue", "Filtered" = "darkgreen")) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "top",
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "gray90"),
      strip.text = ggplot2::element_text(face = "bold")
    )

  if (!is.null(plot_dir)) {
    ggplot2::ggsave(
      filename = file.path(plot_dir, paste0(prefix, "_sparsity_comparison.png")),
      plot = p2,
      width = width,
      height = height * 0.7,
      dpi = dpi
    )
  }

  # === Plot 3: Retention Barplot ===
  total_reads_orig <- sum(orig_abund)
  total_reads_filt <- sum(filt_abund)
  read_retention <- (total_reads_filt / total_reads_orig) * 100
  feature_retention <- (length(filt_features) / length(orig_features)) * 100
  sample_retention <- (ncol(filt_abund) / ncol(orig_abund)) * 100

  retention_data <- data.frame(
    metric = c("Read Retention", "Feature Retention", "Sample Retention"),
    percentage = c(read_retention, feature_retention, sample_retention)
  )

  p3 <- ggplot2::ggplot(retention_data, ggplot2::aes(x = metric, y = percentage, fill = metric)) +
    ggplot2::geom_col(color = "white", linewidth = 0.3) +
    ggplot2::labs(
      title = "Retention Rates After Filtering",
      x = NULL,
      y = "Percentage Retained (%)"
    ) +
    ggplot2::scale_fill_manual(values = rep("steelblue", 3)) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::ylim(0, 100)

  if (!is.null(plot_dir)) {
    ggplot2::ggsave(
      filename = file.path(plot_dir, paste0(prefix, "_retention_rates.png")),
      plot = p3,
      width = width,
      height = height * 0.7,
      dpi = dpi
    )
  }

  # === Plot 5: Heatmap of Top Variable Features ===
  # Calculate relative abundances for heatmap analysis
  orig_rel <- orig_abund / colSums(orig_abund)
  filt_rel <- filt_abund / colSums(filt_abund)

  # Find common features
  common_features <- intersect(orig_features, filt_features)

  p5_original <- NULL
  p5_filtered <- NULL

  if (length(common_features) > 0) {
    # Get top variable features from common set
    orig_common <- orig_rel[match(common_features, orig_features), , drop = FALSE]
    variances <- apply(orig_common, 1, var)
    top_var_idx <- order(variances, decreasing = TRUE)[seq_len(min(heatmap_top, length(variances)))]
    top_var_features <- common_features[top_var_idx]

    # Subset to top variable features and common samples
    common_samples <- intersect(colnames(orig_abund), colnames(filt_abund))
    heatmap_data_orig <- orig_common[top_var_idx, match(common_samples, colnames(orig_common)), drop = FALSE]
    heatmap_data_filt <- filt_rel[match(top_var_features, filt_features),
                                   match(common_samples, colnames(filt_abund)), drop = FALSE]

    # Row-scaled heatmaps
    heatmap_data_orig_scaled <- t(scale(t(heatmap_data_orig)))
    heatmap_data_filt_scaled <- t(scale(t(heatmap_data_filt)))

    # Check if pheatmap is available for better ggplot2 integration
    if (requireNamespace("pheatmap", quietly = TRUE)) {
      # Create descriptive titles
      orig_title <- paste0("Top ", heatmap_top, " Most Variable Features (Original)")
      filt_title <- paste0("Top ", heatmap_top, " Most Variable Features (Filtered)")

      # Calculate appropriate fontsize based on number of features
      row_fontsize <- max(3, min(8, floor(60 / nrow(heatmap_data_orig_scaled))))

      # Original heatmap
      if (!is.null(plot_dir)) {
        pheatmap::pheatmap(
          heatmap_data_orig_scaled,
          main = paste(orig_title, "\n(Rows = Features, Columns = Samples)"),
          color = colorRampPalette(c("blue", "white", "red"))(100),
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          show_rownames = TRUE,
          show_colnames = FALSE,
          fontsize_row = row_fontsize,
          border_color = NA,
          legend = TRUE,
          filename = file.path(plot_dir, paste0(prefix, "_heatmap_original.png")),
          width = width * 1.5,
          height = height * 1.5,
          res = dpi
        )
      }

      # Filtered heatmap
      if (!is.null(plot_dir)) {
        pheatmap::pheatmap(
          heatmap_data_filt_scaled,
          main = paste(filt_title, "\n(Rows = Features, Columns = Samples)"),
          color = colorRampPalette(c("blue", "white", "red"))(100),
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          show_rownames = TRUE,
          show_colnames = FALSE,
          fontsize_row = row_fontsize,
          border_color = NA,
          legend = TRUE,
          filename = file.path(plot_dir, paste0(prefix, "_heatmap_filtered.png")),
          width = width * 1.5,
          height = height * 1.5,
          res = dpi
        )
      }

      # Create ggplot objects for display
      p5_original <- create_heatmap_ggplot(heatmap_data_orig_scaled,
                                            top_var_features,
                                            common_samples,
                                            orig_title)
      p5_filtered <- create_heatmap_ggplot(heatmap_data_filt_scaled,
                                            top_var_features,
                                            common_samples,
                                            filt_title)
    } else {
      # Fallback: try ComplexHeatmap
      if (requireNamespace("ComplexHeatmap", quietly = TRUE) && requireNamespace("grid", quietly = TRUE)) {
        if (!is.null(plot_dir)) {
          # Original heatmap
          ComplexHeatmap::Heatmap(
            heatmap_data_orig_scaled,
            name = paste0(prefix, "_orig"),
            col = colorRampPalette(c("blue", "white", "red"))(100),
            column_title = "Original Table",
            show_row_names = TRUE,
            show_column_names = FALSE,
            row_names_gp = grid::gpar(fontsize = 6),
            column_names_gp = grid::gpar(fontsize = 8)
          )
          png(file.path(plot_dir, paste0(prefix, "_heatmap_original.png")),
              width = width * dpi, height = height * dpi, res = dpi)
          ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
            heatmap_data_orig_scaled,
            name = paste0(prefix, "_orig"),
            col = colorRampPalette(c("blue", "white", "red"))(100),
            column_title = paste("Top", heatmap_top, "Variable Features\nOriginal Table"),
            show_row_names = TRUE,
            show_column_names = FALSE,
            row_names_gp = grid::gpar(fontsize = 6),
            column_names_gp = grid::gpar(fontsize = 8)
          ))
          dev.off()

          # Filtered heatmap
          png(file.path(plot_dir, paste0(prefix, "_heatmap_filtered.png")),
              width = width * dpi, height = height * dpi, res = dpi)
          ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
            heatmap_data_filt_scaled,
            name = paste0(prefix, "_filt"),
            col = colorRampPalette(c("blue", "white", "red"))(100),
            column_title = paste("Top", heatmap_top, "Variable Features\nFiltered Table"),
            show_row_names = TRUE,
            show_column_names = FALSE,
            row_names_gp = grid::gpar(fontsize = 6),
            column_names_gp = grid::gpar(fontsize = 8)
          ))
          dev.off()
        }
      } else {
        warning("Package 'pheatmap' or 'ComplexHeatmap' required for heatmap visualization.")
      }
    }
  }

  # === Plot 7: Top N Features Stacked Barplot ===
  # Calculate total abundance per feature (sum across all samples)
  orig_total_feat <- rowSums(orig_abund)
  filt_total_feat <- rowSums(filt_abund)

  # Convert to relative abundances (within each table)
  orig_rel_total <- orig_total_feat / sum(orig_total_feat)
  filt_rel_total <- filt_total_feat / sum(filt_total_feat)

  # Get top N features independently from each table
  orig_top_idx <- order(orig_rel_total, decreasing = TRUE)[seq_len(min(top_n, length(orig_rel_total)))]
  filt_top_idx <- order(filt_rel_total, decreasing = TRUE)[seq_len(min(top_n, length(filt_rel_total)))]

  orig_top_features <- orig_features[orig_top_idx]
  filt_top_features <- filt_features[filt_top_idx]

  orig_top_rels <- orig_rel_total[orig_top_idx]
  filt_top_rels <- filt_rel_total[filt_top_idx]

  # Create combined feature list (union of top N from both tables)
  all_top_features <- unique(c(orig_top_features, filt_top_features))
  n_plot_features <- length(all_top_features)

  # Build relative abundance vectors for plotting
  orig_plot_rels <- vapply(all_top_features, function(f) {
    if (f %in% orig_top_features) {
      orig_top_rels[match(f, orig_top_features)]
    } else {
      0
    }
  }, numeric(1))

  filt_plot_rels <- vapply(all_top_features, function(f) {
    if (f %in% filt_top_features) {
      filt_top_rels[match(f, filt_top_features)]
    } else {
      0
    }
  }, numeric(1))

  # Calculate "Other" category
  orig_other <- 1 - sum(orig_plot_rels)
  filt_other <- 1 - sum(filt_plot_rels)

  plot_features <- c(all_top_features, "Other")
  n_total_segments <- length(plot_features)  # includes "Other"

  # Generate qualitative colors for all segments (including Other)
  qual_colors <- c(
    "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
    "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#E5C494"
  )

  if (n_total_segments <= length(qual_colors)) {
    plot_colors <- qual_colors[seq_len(n_total_segments)]
  } else {
    # Extend by cycling through the palette
    plot_colors <- qual_colors[((seq_len(n_total_segments) - 1) %% length(qual_colors)) + 1]
  }

  # Last color is for "Other" (gray)
  plot_colors[n_total_segments] <- "#999999"

  # Build final vectors including "Other"
  orig_final_rels <- c(orig_plot_rels, orig_other)
  filt_final_rels <- c(filt_plot_rels, filt_other)

  # Create data frame for stacked barplot
  stacked_data <- data.frame(
    feature = factor(rep(plot_features, 2), levels = rev(plot_features)),
    table_type = rep(c("Original", "Filtered"), each = length(plot_features)),
    proportion = c(orig_final_rels, filt_final_rels)
  )

  # Create named color vector
  fill_colors <- setNames(plot_colors, plot_features)

  p6 <- ggplot2::ggplot(stacked_data, ggplot2::aes(x = table_type, y = proportion, fill = feature)) +
    ggplot2::geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = fill_colors, na.value = "#999999") +
    ggplot2::labs(
      title = sprintf("Top %d Features by Total Relative Abundance", top_n),
      x = NULL,
      y = "Relative Abundance",
      fill = "Feature"
    ) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      panel.grid.minor.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 10)
    )

  if (!is.null(plot_dir)) {
    ggplot2::ggsave(
      filename = file.path(plot_dir, paste0(prefix, "_top_features_stacked.png")),
      plot = p6,
      width = width,
      height = height,
      dpi = dpi
    )
  }

  # Return info about saved plots and ggplot objects
  result <- list(
    plots = list(
      sample_coverage_distribution = p1,
      feature_abundance_distribution = p_feat,
      sparsity_comparison = p2,
      retention_rates = p3,
      top_features_stacked = p6,
      heatmap_original = p5_original,
      heatmap_filtered = p5_filtered
    )
  )

  if (!is.null(plot_dir)) {
    result$paths <- list(
      sample_coverage_distribution = file.path(plot_dir, paste0(prefix, "_sample_coverage_distribution.png")),
      feature_abundance_distribution = file.path(plot_dir, paste0(prefix, "_feature_abundance_distribution.png")),
      sparsity_comparison = file.path(plot_dir, paste0(prefix, "_sparsity_comparison.png")),
      retention_rates = file.path(plot_dir, paste0(prefix, "_retention_rates.png")),
      top_features_stacked = file.path(plot_dir, paste0(prefix, "_top_features_stacked.png")),
      heatmap_original = file.path(plot_dir, paste0(prefix, "_heatmap_original.png")),
      heatmap_filtered = file.path(plot_dir, paste0(prefix, "_heatmap_filtered.png"))
    )
    return(result)
  }

  invisible(result)
}

# Helper function to create ggplot heatmap from matrix
create_heatmap_ggplot <- function(mat, row_names, col_names, main_title) {
  # Convert matrix to long format
  mat_df <- as.data.frame(mat)
  mat_df$feature <- row_names
  mat_long <- tidyr::pivot_longer(mat_df,
                                   cols = -feature,
                                   names_to = "sample",
                                   values_to = "value")

  # Add sample index for ordering
  mat_long$sample_idx <- match(mat_long$sample, col_names)

  # Calculate number of sample ticks to show (more for smaller datasets)
  n_samples <- length(col_names)
  if (n_samples <= 10) {
    sample_tick_interval <- 1
  } else if (n_samples <= 50) {
    sample_tick_interval <- max(1, floor(n_samples / 8))
  } else {
    sample_tick_interval <- max(1, floor(n_samples / 10))
  }

  # Factor levels to maintain order
  mat_long$feature <- factor(mat_long$feature, levels = rev(unique(mat_long$feature)))

  # Create ggplot heatmap
  p <- ggplot2::ggplot(mat_long, ggplot2::aes(x = sample_idx, y = feature, fill = value)) +
    ggplot2::geom_tile(color = "gray90", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      name = "Z-score\n(scaled)",
      limits = c(-2, 2),
      oob = scales::squish
    ) +
    ggplot2::labs(
      title = main_title,
      x = "Samples",
      y = "Features"
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(1, n_samples, by = sample_tick_interval),
      labels = col_names[seq(1, n_samples, by = sample_tick_interval)]
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      axis.text.y = ggplot2::element_text(size = 8, family = "monospace"),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right",
      legend.box = "vertical",
      aspect.ratio = 0.6
    )

  return(p)
}
