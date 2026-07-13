# =============================================================================
# Pipeline Reporting Helper Functions (Internal)
# =============================================================================
# These functions generate QC plots and text reports for the pipeline.
# All functions are internal (. prefix).
# =============================================================================


#' Generate all QC comparison plots
#'
#' Creates and saves all QC visualization plots comparing original and filtered tables.
#'
#' @param original_table Original feature table
#' @param filtered_table Filtered feature table
#' @param output_dir Output directory for plots
#' @param prefix Prefix for output filenames
#' @param verbose Logical. Print progress messages?
#'
#' @return List of paths to generated plot files
.generate_qc_plots <- function(original_table, filtered_table, output_dir, prefix, verbose = TRUE) {
  plot_paths <- list()

  if (verbose) cat("Generating QC comparison plots...\n")

  # Main QC comparison plots
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plot_paths$qc_comparison <- plot_qc_comparison(
      original_table, filtered_table,
      plot_dir = output_dir, prefix = prefix
    )
  }

  plot_paths
}


#' Generate sparsity elbow plot
#'
#' Saves a plot of the sparsity elbow detection result.
#'
#' @param elbow_result Result from identify_sparsity_elbow()
#' @param output_dir Output directory
#' @param prefix Prefix for filename
#' @param main Plot title
#' @param verbose Logical. Print progress messages?
#'
#' @return Path to saved plot file
.generate_sparsity_elbow_plot <- function(elbow_result, output_dir, prefix, main, verbose = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 required for sparsity elbow plot")
    return(NULL)
  }

  if (verbose) cat("Generating Sparsity Elbow plot...\n")
  p <- plot_sparsity_elbow(elbow_result, main = main)
  plot_path <- file.path(output_dir, paste0(prefix, "_sparsity_elbow.png"))
  ggplot2::ggsave(plot_path, plot = p, width = 10, height = 8, dpi = 300)
  if (verbose) cat(sprintf("Sparsity Elbow plot saved to: %s\n", plot_path))
  plot_path
}


#' Generate depth-sparsity outlier plot
#'
#' Saves a plot of the depth-sparsity outlier analysis result.
#'
#' @param analysis_result Result from analyze_depth_sparsity()
#' @param output_dir Output directory
#' @param prefix Prefix for filename
#' @param main Plot title
#' @param verbose Logical. Print progress messages?
#'
#' @return Path to saved plot file
.generate_depth_sparsity_plot <- function(analysis_result, output_dir, prefix, main, verbose = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 required for depth-sparsity plot")
    return(NULL)
  }

  if (verbose) cat("Generating Depth-Sparsity Outlier plot...\n")
  p <- plot_depth_sparsity(analysis_result, main = main)
  plot_path <- file.path(output_dir, paste0(prefix, "_depth_sparsity_outliers.png"))
  ggplot2::ggsave(plot_path, plot = p, width = 10, height = 8, dpi = 300)
  if (verbose) cat(sprintf("Depth-Sparsity plot saved to: %s\n", plot_path))
  plot_path
}


#' Generate scree analysis plot
#'
#' Saves a plot of the scree/saturation analysis result.
#'
#' @param scree_result Result from compute_scree()
#' @param output_dir Output directory
#' @param prefix Prefix for filename
#' @param main Plot title
#' @param verbose Logical. Print progress messages?
#'
#' @return Path to saved plot file
.generate_scree_plot <- function(scree_result, output_dir, prefix, main, verbose = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 required for scree plot")
    return(NULL)
  }

  if (verbose) cat("Generating Scree/Saturation plot...\n")
  p <- plot_scree(scree_result, main = main)
  plot_path <- file.path(output_dir, paste0(prefix, "_scree_analysis.png"))
  ggplot2::ggsave(plot_path, plot = p, width = 10, height = 8, dpi = 300)
  if (verbose) cat(sprintf("Scree plot saved to: %s\n", plot_path))
  plot_path
}


#' Generate text summary report
#'
#' Creates a comprehensive text report of the filtering pipeline results.
#'
#' @param original_stats Stats before filtering
#' @param qc_metrics QC metrics from compute_filtering_qc()
#' @param presence_stats Presence analysis statistics
#' @param filtering_steps List of filtering step summaries
#' @param sparsity_elbow_result Sparsity elbow result (optional)
#' @param depth_sparsity_result Depth-sparsity result (optional)
#' @param scree_result Scree analysis result (optional)
#' @param output_dir Output directory
#' @param prefix Prefix for filename
#' @param verbose Logical. Print progress messages?
#'
#' @return Path to saved report file
.generate_filtering_report <- function(original_stats, qc_metrics, presence_stats,
                                        filtering_steps, sparsity_elbow_result,
                                        depth_sparsity_result, scree_result,
                                        output_dir, prefix, verbose = TRUE) {
  if (verbose) cat("Generating summary report...\n")

  report_lines <- c(
    "========================================",
    "FEATURE TABLE FILTERING PIPELINE REPORT",
    "========================================",
    "",
    "--- ORIGINAL TABLE STATISTICS ---",
    sprintf("Features: %d", original_stats$features),
    sprintf("Samples: %d", original_stats$samples),
    sprintf("Total Reads: %d", original_stats$reads),
    ""
  )

  # Add filtering steps summary
  if (length(filtering_steps) > 0) {
    report_lines <- c(report_lines,
      "--- FILTERING STEPS ---"
    )
    for (step_name in names(filtering_steps)) {
      step <- filtering_steps[[step_name]]
      report_lines <- c(report_lines,
        sprintf("\n%s:", toupper(step_name)),
        sprintf("  Method: %s", step$method),
        sprintf("  Samples: %d -> %d (removed: %d)",
                step$samples_before, step$samples_after, step$samples_removed),
        sprintf("  Features: %d -> %d (removed: %d)",
                step$features_before, step$features_after, step$features_removed),
        sprintf("  Reads: %d -> %d (removed: %d)",
                step$reads_before, step$reads_after, step$reads_removed)
      )
    }
    report_lines <- c(report_lines, "")
  }

  # Add QC metrics summary
  report_lines <- c(report_lines,
    "--- QC METRICS SUMMARY ---",
    sprintf("Sparsity: %.2f%% -> %.2f%% (drop: %.2f pp)",
            qc_metrics$sparsity_original * 100,
            qc_metrics$sparsity_filtered * 100,
            qc_metrics$sparsity_drop_percent),
    sprintf("Read Retention: %.2f%%", qc_metrics$read_retention_percent),
    sprintf("Feature Retention: %.2f%%", qc_metrics$feature_retention_percent),
    sprintf("Sample Retention: %.2f%%", qc_metrics$sample_retention_percent),
    ""
  )

  # Add Shannon ENS retention
  if (!is.na(qc_metrics$shannon_ens_retention_percent)) {
    report_lines <- c(report_lines,
      sprintf("Shannon ENS Retention: %.2f%%", qc_metrics$shannon_ens_retention_percent)
    )
  }

  # Add Simpson ENS retention
  if (!is.na(qc_metrics$simpson_ens_retention_percent)) {
    report_lines <- c(report_lines,
      sprintf("Simpson ENS Retention: %.2f%%", qc_metrics$simpson_ens_retention_percent)
    )
  }

  report_lines <- c(report_lines,
    "",
    "--- TOP N FEATURE COMPARISON ---",
    sprintf("Top %d Overlap: %d features (%.2f%%)",
            length(qc_metrics$orig_top_features),
            qc_metrics$top_n_overlap_count,
            qc_metrics$top_n_overlap_percent),
    sprintf("Jaccard Similarity: %.4f", qc_metrics$top_n_jaccard_similarity),
    sprintf("Rank Correlation: %.4f (p = %.4f)",
            qc_metrics$rank_abundance_correlation,
            qc_metrics$rank_abundance_pvalue),
    "",
    "========================================",
    "END OF REPORT",
    "========================================"
  )

  # Save report
  report_path <- file.path(output_dir, paste0(prefix, "_filtering_report.txt"))
  writeLines(report_lines, report_path)
  if (verbose) cat(sprintf("Report saved to: %s\n", report_path))
  report_path
}
