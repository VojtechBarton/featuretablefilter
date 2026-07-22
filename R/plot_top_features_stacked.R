#' Plot top N features as stacked barplots comparing original and filtered tables using ggplot2
#'
#' Creates a visualization showing the relative abundance composition of the top N features
#' in both original and filtered feature tables. Each table's top N features are identified
#' independently by summing abundances across all samples, converting to relative abundances,
#' and selecting the top N. Non-top-N features are grouped as "Other".
#'
#' @param original_table Original feature table (data.frame) before filtering.
#' @param filtered_table Filtered feature table (data.frame) after filtering.
#' @param top_n Number of top features to show. Default is 10.
#' @param plot_dir Optional directory path to save the plot. If NULL, plot is displayed only.
#' @param prefix Prefix for saved plot filename. Default is "top_features_stacked".
#' @param colors Optional vector of colors for features. If NULL, auto-generated using ColorBrewer-style palette.
#' @param width Plot width in inches for ggsave. Default 10.
#' @param height Plot height in inches for ggsave. Default 7.
#' @param dpi DPI for saved raster output. Default 300.
#' @param main_title Main title for the plot. Default is auto-generated.
#'
#' @return Returns a list containing:
#' \describe{
#'   \item{plot_path}{Path to saved plot (if plot_dir provided)}
#'   \item{orig_top_features}{Character vector of top N features in original table}
#'   \item{filt_top_features}{Character vector of top N features in filtered table}
#'   \item{orig_rels}{Named vector of relative abundances for original top features + Other}
#'   \item{filt_rels}{Named vector of relative abundances for filtered top features + Other}
#'   \item{plot}{The ggplot object}
#' }
#'
#' @export
#'
#' @examples
#' data(example_feature_table)
#' filtered <- filter_by_coverage(example_feature_table, min_reads = 1000)
#' plot_top_features_stacked(example_feature_table, filtered, top_n = 5)
plot_top_features_stacked <- function(original_table, filtered_table,
                                       top_n = 10, plot_dir = NULL,
                                       prefix = "top_features_stacked",
                                       colors = NULL,
                                       width = 10, height = 7, dpi = 300,
                                       main_title = NULL) {
  # Extract abundance matrices
  orig_abund <- as.matrix(original_table[, -1, drop = FALSE])
  filt_abund <- as.matrix(filtered_table[, -1, drop = FALSE])
  orig_features <- original_table[, 1]
  filt_features <- filtered_table[, 1]

  # Calculate total abundance per feature (sum across all samples)
  orig_total <- rowSums(orig_abund)
  filt_total <- rowSums(filt_abund)

  # Convert to relative abundances (within each table)
  orig_rel_total <- orig_total / sum(orig_total)
  filt_rel_total <- filt_total / sum(filt_total)

  # Get top N features independently from each table
  orig_top_idx <- order(orig_rel_total, decreasing = TRUE)[seq_len(min(top_n, length(orig_rel_total)))]
  filt_top_idx <- order(filt_rel_total, decreasing = TRUE)[seq_len(min(top_n, length(filt_rel_total)))]

  orig_top_features <- orig_features[orig_top_idx]
  filt_top_features <- filt_features[filt_top_idx]

  orig_top_rels <- orig_rel_total[orig_top_idx]
  filt_top_rels <- filt_rel_total[filt_top_idx]

  # Create combined feature list (union of top N from both tables)
  all_top_features <- unique(c(orig_top_features, filt_top_features))
  n_features <- length(all_top_features)

  # Build relative abundance vectors for plotting (aligned to all_top_features)
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

  # Calculate "Other" category (remaining abundance)
  orig_other <- 1 - sum(orig_plot_rels)
  filt_other <- 1 - sum(filt_plot_rels)

  # Add "Other" to the feature list for plotting
  plot_features <- c(all_top_features, "Other")
  n_total_segments <- length(plot_features)
  orig_final_rels <- c(orig_plot_rels, orig_other)
  filt_final_rels <- c(filt_plot_rels, filt_other)

  # Generate qualitative colors if not provided
  if (is.null(colors)) {
    # Use ColorBrewer Set3 qualitative palette (12 distinct colors) + extended
    qual_colors <- c(
      "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
      "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
      "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#E5C494"
    )

    if (n_total_segments <= length(qual_colors)) {
      colors <- qual_colors[seq_len(n_total_segments)]
    } else {
      # Extend by cycling through the palette
      colors <- qual_colors[((seq_len(n_total_segments) - 1) %% length(qual_colors)) + 1]
    }
  }

  # Last color is for "Other" (gray)
  plot_colors <- colors
  plot_colors[n_total_segments] <- "#999999"

  # Create data frame for ggplot (long format for stacked bars)
  plot_data <- data.frame(
    feature = rep(plot_features, 2),
    table_type = rep(c("Original", "Filtered"), each = length(plot_features)),
    proportion = c(orig_final_rels, filt_final_rels)
  )

  # Create named color vector for scale_fill_manual
  fill_colors <- setNames(plot_colors, plot_features)

  # Set default title
  if (is.null(main_title)) {
    main_title <- sprintf("Top %d Features by Total Relative Abundance", top_n)
  }

  # Build ggplot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = table_type, y = proportion, fill = feature)) +
    ggplot2::geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = fill_colors, na.value = "#999999") +
    ggplot2::labs(
      title = main_title,
      x = NULL,
      y = "Relative Abundance",
      fill = "Feature"
    ) +
    ggplot2::coord_flip() +  # Horizontal bars for better label readability
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      panel.grid.minor.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 10)
    )

  # Add annotation about independent top N selection
  annotation_text <- sprintf("* Top %d features selected independently for each table\n* Based on sum of abundances across all samples converted to relative abundances", top_n)

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
    plot_path <- file.path(plot_dir, paste0(prefix, ".png"))
  } else {
    plot_path <- NULL
  }

  return(list(
    plot_path = plot_path,
    orig_top_features = orig_top_features,
    filt_top_features = filt_top_features,
    orig_rels = setNames(orig_final_rels, plot_features),
    filt_rels = setNames(filt_final_rels, plot_features),
    plot = p
  ))
}
