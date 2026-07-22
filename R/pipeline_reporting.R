# =============================================================================
# Pipeline Reporting Helper Functions (Internal)
# =============================================================================
# These functions generate QC plots and text reports for the pipeline.
# All functions are internal (. prefix).
# =============================================================================

#' Generate all QC comparison plots#'#' Creates and saves all QC visualization plots comparing original and filtered tables.#'#' @param original_table Original feature table#' @param filtered_table Filtered feature table#' @param output_dir Output directory for plots#' @param prefix Prefix for output filenames#' @param verbose Logical. Print progress messages?#'#' @return List of paths to generated plot files
#' @noRd
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

#' Generate sparsity elbow plot#'#' Saves a plot of the sparsity elbow detection result.#'#' @param elbow_result Result from identify_sparsity_elbow()#' @param output_dir Output directory#' @param prefix Prefix for filename#' @param main Plot title#' @param verbose Logical. Print progress messages?#'#' @return Path to saved plot file
#' @noRd
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

#' Generate depth-sparsity outlier plot#'#' Saves a plot of the depth-sparsity outlier analysis result.#'#' @param analysis_result Result from analyze_depth_sparsity()#' @param output_dir Output directory#' @param prefix Prefix for filename#' @param main Plot title#' @param verbose Logical. Print progress messages?#'#' @return Path to saved plot file
#' @noRd
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

#' Generate scree analysis plot#'#' Saves a plot of the scree/saturation analysis result.#'#' @param scree_result Result from compute_scree()#' @param output_dir Output directory#' @param prefix Prefix for filename#' @param main Plot title#' @param verbose Logical. Print progress messages?#'#' @return Path to saved plot file
#' @noRd
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

#' Generate text summary report#'#' Creates a comprehensive text report of the filtering pipeline results.#'#' @param original_stats Stats before filtering#' @param qc_metrics QC metrics from compute_filtering_qc()#' @param presence_stats Presence analysis statistics#' @param filtering_steps List of filtering step summaries#' @param sparsity_elbow_result Sparsity elbow result (optional)#' @param depth_sparsity_result Depth-sparsity result (optional)#' @param scree_result Scree analysis result (optional)#' @param output_dir Output directory#' @param prefix Prefix for filename#' @param verbose Logical. Print progress messages?#' @param input_description Description of input source (file path or object type)#' @param pipeline_params Named list of pipeline parameters used#' @param filtered_stats Stats after filtering (optional, for comparison)#'#' @return Path to saved report file
#' @noRd
.generate_filtering_report <- function(original_stats, qc_metrics, presence_stats,
                                        filtering_steps, sparsity_elbow_result,
                                        depth_sparsity_result, scree_result,
                                        output_dir, prefix, verbose = TRUE,
                                        input_description = NULL,
                                        pipeline_params = NULL,
                                        filtered_stats = NULL) {
  if (verbose) cat("Generating summary report...\n")

  # Generate header with metadata
  report_lines <- c(
    "================================================================================",
    "                    FEATURE TABLE FILTERING PIPELINE REPORT",
    "================================================================================",
    "",
    sprintf("Report Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("Pipeline Version: %s", utils::packageVersion("featuretablefilter")),
    ""
  )

  # Input and output information
  report_lines <- c(report_lines,
    "--- INPUT / OUTPUT ---",
    sprintf("Input Source: %s", input_description %||% "Unknown"),
    sprintf("Output Prefix: %s", prefix),
    sprintf("Output Directory: %s", normalizePath(output_dir, mustWork = FALSE)),
    ""
  )

  # Pipeline parameters summary
  if (!is.null(pipeline_params) && length(pipeline_params) > 0) {
    report_lines <- c(report_lines, "--- PIPELINE PARAMETERS ---")
    param_lines <- vapply(names(pipeline_params), function(name) {
      val <- pipeline_params[[name]]
      if (is.null(val)) {
        sprintf("  %s: NULL", name)
      } else if (is.logical(val)) {
        sprintf("  %s: %s", name, ifelse(val, "TRUE", "FALSE"))
      } else if (is.numeric(val) && length(val) == 1) {
        sprintf("  %s: %.4g", name, val)
      } else if (is.character(val) && length(val) == 1) {
        sprintf("  %s: \"%s\"", name, val)
      } else {
        sprintf("  %s: %s", name, paste(val, collapse = ", "))
      }
    }, character(1))
    report_lines <- c(report_lines, param_lines, "")
  }

  report_lines <- c(report_lines,
    "--- ORIGINAL TABLE STATISTICS ---",
    sprintf("Features: %d", original_stats$features),
    sprintf("Samples: %d", original_stats$samples),
    sprintf("Total Reads: %d", original_stats$reads),
    ""
  )

  # Add Hill numbers (diversity) for original table
  if (!is.na(qc_metrics$shannon_ens_original)) {
    report_lines <- c(report_lines,
      sprintf("Shannon ENS (q=1): %.2f", qc_metrics$shannon_ens_original)
    )
  }
  if (!is.na(qc_metrics$simpson_ens_original)) {
    report_lines <- c(report_lines,
      sprintf("Simpson ENS (q=2): %.2f", qc_metrics$simpson_ens_original)
    )
  }
  report_lines <- c(report_lines, "")

  # Add filtered table statistics if available
  if (!is.null(filtered_stats)) {
    report_lines <- c(report_lines,
      "--- FILTERED TABLE STATISTICS ---",
      sprintf("Features: %d (removed: %d)",
              filtered_stats$features, original_stats$features - filtered_stats$features),
      sprintf("Samples: %d (removed: %d)",
              filtered_stats$samples, original_stats$samples - filtered_stats$samples),
      sprintf("Total Reads: %d (removed: %d)",
              filtered_stats$reads, original_stats$reads - filtered_stats$reads),
      ""
    )

    # Add Hill numbers for filtered table
    if (!is.na(qc_metrics$shannon_ens_filtered)) {
      report_lines <- c(report_lines,
        sprintf("Shannon ENS (q=1): %.2f (retained: %.2f%%)",
                qc_metrics$shannon_ens_filtered, qc_metrics$shannon_ens_retention_percent)
      )
    }
    if (!is.na(qc_metrics$simpson_ens_filtered)) {
      report_lines <- c(report_lines,
        sprintf("Simpson ENS (q=2): %.2f (retained: %.2f%%)",
                qc_metrics$simpson_ens_filtered, qc_metrics$simpson_ens_retention_percent)
      )
    }
    report_lines <- c(report_lines, "")
  }

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

  # Add Procrustes analysis (compositional similarity)
  if (!is.na(qc_metrics$procrustes_m2)) {
    report_lines <- c(report_lines,
      ""
    )
    if (!is.na(qc_metrics$procrustes_correlation)) {
      report_lines <- c(report_lines,
        sprintf("Procrustes Correlation: %.4f", qc_metrics$procrustes_correlation)
      )
      if (!is.na(qc_metrics$procrustes_pvalue)) {
        report_lines <- c(report_lines,
          sprintf("Procrustes p-value: %.4f", qc_metrics$procrustes_pvalue)
        )
      }
    }
    report_lines <- c(report_lines,
      sprintf("Procrustes M^2 Statistic: %.6f", qc_metrics$procrustes_m2),
      "(Lower M^2 = more similar compositional structure between original and filtered)"
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

  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save text report
  report_path <- file.path(output_dir, paste0(prefix, "_filtering_report.txt"))
  writeLines(report_lines, report_path)
  if (verbose) cat(sprintf("Text report saved to: %s\n", report_path))

  # Generate Markdown report
  md_path <- .generate_markdown_report(
    original_stats, qc_metrics, presence_stats, filtering_steps,
    sparsity_elbow_result, depth_sparsity_result, scree_result,
    output_dir, prefix, input_description, pipeline_params, filtered_stats
  )
  if (!is.null(md_path) && verbose) {
    cat(sprintf("Markdown report saved to: %s\n", md_path))
  }

  # Generate PDF report
  pdf_path <- .generate_pdf_report(
    original_stats, qc_metrics, presence_stats, filtering_steps,
    sparsity_elbow_result, depth_sparsity_result, scree_result,
    output_dir, prefix, input_description, pipeline_params, filtered_stats, verbose
  )
  if (!is.null(pdf_path) && verbose) {
    cat(sprintf("PDF report saved to: %s\n", pdf_path))
  }

  list(text = report_path, markdown = md_path, pdf = pdf_path)
}

#' Generate Markdown version of the filtering report
#' @noRd
.generate_markdown_report <- function(original_stats, qc_metrics, presence_stats,
                                       filtering_steps, sparsity_elbow_result,
                                       depth_sparsity_result, scree_result,
                                       output_dir, prefix, input_description,
                                       pipeline_params, filtered_stats) {
  md_lines <- c(
    "# Feature Table Filtering Pipeline Report",
    "",
    sprintf("**Report Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("**Pipeline Version:** %s", utils::packageVersion("featuretablefilter")),
    ""
  )

  # Input and output
  md_lines <- c(md_lines,
    "## Input / Output",
    "",
    sprintf("- **Input Source:** %s", input_description %||% "Unknown"),
    sprintf("- **Output Prefix:** %s", prefix),
    sprintf("- **Output Directory:** %s", normalizePath(output_dir, mustWork = FALSE)),
    ""
  )

  # Pipeline parameters
  if (!is.null(pipeline_params) && length(pipeline_params) > 0) {
    md_lines <- c(md_lines, "## Pipeline Parameters", "")
    md_lines <- c(md_lines, "| Parameter | Value |")
    md_lines <- c(md_lines, "|-----------|-------|")
    for (name in names(pipeline_params)) {
      val <- pipeline_params[[name]]
      if (is.null(val)) {
        val_str <- "NULL"
      } else if (is.logical(val)) {
        val_str <- ifelse(val, "TRUE", "FALSE")
      } else if (is.numeric(val) && length(val) == 1) {
        val_str <- sprintf("%.4g", val)
      } else if (is.character(val) && length(val) == 1) {
        val_str <- sprintf("`%s`", val)
      } else {
        val_str <- paste(val, collapse = ", ")
      }
      md_lines <- c(md_lines, sprintf("| %s | %s |", name, val_str))
    }
    md_lines <- c(md_lines, "")
  }

  # Original table statistics
  md_lines <- c(md_lines,
    "## Original Table Statistics",
    "",
    sprintf("- **Features:** %d", original_stats$features),
    sprintf("- **Samples:** %d", original_stats$samples),
    sprintf("- **Total Reads:** %d", original_stats$reads),
    ""
  )

  # Hill numbers for original
  if (!is.na(qc_metrics$shannon_ens_original)) {
    md_lines <- c(md_lines, sprintf("- **Shannon ENS (q=1):** %.2f", qc_metrics$shannon_ens_original))
  }
  if (!is.na(qc_metrics$simpson_ens_original)) {
    md_lines <- c(md_lines, sprintf("- **Simpson ENS (q=2):** %.2f", qc_metrics$simpson_ens_original))
  }
  md_lines <- c(md_lines, "")

  # Filtered table statistics
  if (!is.null(filtered_stats)) {
    md_lines <- c(md_lines,
      "## Filtered Table Statistics",
      "",
      sprintf("- **Features:** %d (*removed:* %d)",
              filtered_stats$features, original_stats$features - filtered_stats$features),
      sprintf("- **Samples:** %d (*removed:* %d)",
              filtered_stats$samples, original_stats$samples - filtered_stats$samples),
      sprintf("- **Total Reads:** %d (*removed:* %d)",
              filtered_stats$reads, original_stats$reads - filtered_stats$reads),
      ""
    )

    # Hill numbers for filtered
    if (!is.na(qc_metrics$shannon_ens_filtered)) {
      md_lines <- c(md_lines,
        sprintf("- **Shannon ENS (q=1):** %.2f (*retained:* %.2f%%)",
                qc_metrics$shannon_ens_filtered, qc_metrics$shannon_ens_retention_percent)
      )
    }
    if (!is.na(qc_metrics$simpson_ens_filtered)) {
      md_lines <- c(md_lines,
        sprintf("- **Simpson ENS (q=2):** %.2f (*retained:* %.2f%%)",
                qc_metrics$simpson_ens_filtered, qc_metrics$simpson_ens_retention_percent)
      )
    }
    md_lines <- c(md_lines, "")
  }

  # Filtering steps
  if (length(filtering_steps) > 0) {
    md_lines <- c(md_lines, "## Filtering Steps", "")
    for (step_name in names(filtering_steps)) {
      step <- filtering_steps[[step_name]]
      md_lines <- c(md_lines,
        sprintf("### %s (%s)", toupper(step_name), step$method),
        "",
        sprintf("| Metric | Before | After | Removed |"),
        sprintf("|--------|--------|-------|---------|"),
        sprintf("| Samples | %d | %d | %d |",
                step$samples_before, step$samples_after, step$samples_removed),
        sprintf("| Features | %d | %d | %d |",
                step$features_before, step$features_after, step$features_removed),
        sprintf("| Reads | %d | %d | %d |",
                step$reads_before, step$reads_after, step$reads_removed),
        ""
      )
    }
  }

  # QC metrics summary
  md_lines <- c(md_lines,
    "## QC Metrics Summary",
    "",
    sprintf("- **Sparsity:** %.2f%% -> %.2f%% (drop: %.2f pp)",
            qc_metrics$sparsity_original * 100,
            qc_metrics$sparsity_filtered * 100,
            qc_metrics$sparsity_drop_percent),
    sprintf("- **Read Retention:** %.2f%%", qc_metrics$read_retention_percent),
    sprintf("- **Feature Retention:** %.2f%%", qc_metrics$feature_retention_percent),
    sprintf("- **Sample Retention:** %.2f%%", qc_metrics$sample_retention_percent),
    ""
  )

  # Diversity retention
  if (!is.na(qc_metrics$shannon_ens_retention_percent)) {
    md_lines <- c(md_lines, sprintf("- **Shannon ENS Retention:** %.2f%%", qc_metrics$shannon_ens_retention_percent))
  }
  if (!is.na(qc_metrics$simpson_ens_retention_percent)) {
    md_lines <- c(md_lines, sprintf("- **Simpson ENS Retention:** %.2f%%", qc_metrics$simpson_ens_retention_percent))
  }

  # Procrustes
  if (!is.na(qc_metrics$procrustes_m2)) {
    md_lines <- c(md_lines, "")
    if (!is.na(qc_metrics$procrustes_correlation)) {
      md_lines <- c(md_lines, sprintf("- **Procrustes Correlation:** %.4f", qc_metrics$procrustes_correlation))
    }
    if (!is.na(qc_metrics$procrustes_pvalue)) {
      md_lines <- c(md_lines, sprintf("- **Procrustes p-value:** %.4f", qc_metrics$procrustes_pvalue))
    }
    md_lines <- c(md_lines,
      sprintf("- **Procrustes M^2 Statistic:** %.6f *(lower = more similar)*", qc_metrics$procrustes_m2)
    )
  }
  md_lines <- c(md_lines, "")

  # Top N comparison
  md_lines <- c(md_lines,
    "## Top N Feature Comparison",
    "",
    sprintf("- **Top %d Overlap:** %d features (%.2f%%)",
            length(qc_metrics$orig_top_features),
            qc_metrics$top_n_overlap_count,
            qc_metrics$top_n_overlap_percent),
    sprintf("- **Jaccard Similarity:** %.4f", qc_metrics$top_n_jaccard_similarity),
    sprintf("- **Rank Correlation:** %.4f (p = %.4f)",
            qc_metrics$rank_abundance_correlation,
            qc_metrics$rank_abundance_pvalue),
    ""
  )

  # Save Markdown
  md_path <- file.path(output_dir, paste0(prefix, "_filtering_report.md"))
  writeLines(md_lines, md_path)
  md_path
}

#' Generate PDF version of the filtering report
#' @noRd
.generate_pdf_report <- function(original_stats, qc_metrics, presence_stats,
                                  filtering_steps, sparsity_elbow_result,
                                  depth_sparsity_result, scree_result,
                                  output_dir, prefix, input_description,
                                  pipeline_params, filtered_stats, verbose = TRUE) {
  # Check for required packages
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    if (verbose) cat("  rmarkdown not available - skipping PDF report\n")
    return(NULL)
  }

  # Create temporary Rmd file
  tmp_rmd <- tempfile(fileext = ".Rmd")

  # Generate Rmd content
  rmd_content <- c(
    "---",
    sprintf("title: \"Feature Table Filtering Report\""),
    sprintf("subtitle: \"%s\"", Sys.time()),
    "output:",
    "  pdf_document:",
    "    toc: true",
    "    number_sections: true",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE)",
    "```",
    ""
  )

  # Header
  rmd_content <- c(rmd_content,
    "# Summary",
    "",
    sprintf("- **Input Source:** %s", input_description %||% "Unknown"),
    sprintf("- **Output Prefix:** %s", prefix),
    sprintf("- **Pipeline Version:** %s", utils::packageVersion("featuretablefilter")),
    ""
  )

  # Original table
  rmd_content <- c(rmd_content,
    "## Original Table",
    "",
    sprintf("| Metric | Value |"),
    "|--------|-------|",
    sprintf("| Features | %d |", original_stats$features),
    sprintf("| Samples | %d |", original_stats$samples),
    sprintf("| Total Reads | %d |", original_stats$reads),
    ""
  )

  # Filtered table
  if (!is.null(filtered_stats)) {
    rmd_content <- c(rmd_content,
      "## Filtered Table",
      "",
      sprintf("| Metric | Value | Removed |"),
      "|--------|-------|---------|",
      sprintf("| Features | %d | %d |",
              filtered_stats$features, original_stats$features - filtered_stats$features),
      sprintf("| Samples | %d | %d |",
              filtered_stats$samples, original_stats$samples - filtered_stats$samples),
      sprintf("| Total Reads | %d | %d |",
              filtered_stats$reads, original_stats$reads - filtered_stats$reads),
      ""
    )
  }

  # QC metrics
  rmd_content <- c(rmd_content,
    "## QC Metrics",
    "",
    sprintf("| Metric | Value |"),
    "|--------|-------|",
    sprintf("| Sparsity (original) | %.2f%% |", qc_metrics$sparsity_original * 100),
    sprintf("| Sparsity (filtered) | %.2f%% |", qc_metrics$sparsity_filtered * 100),
    sprintf("| Read Retention | %.2f%% |", qc_metrics$read_retention_percent),
    sprintf("| Feature Retention | %.2f%% |", qc_metrics$feature_retention_percent),
    sprintf("| Sample Retention | %.2f%% |", qc_metrics$sample_retention_percent),
    ""
  )

  # Diversity
  if (!is.na(qc_metrics$shannon_ens_retention_percent) || !is.na(qc_metrics$simpson_ens_retention_percent)) {
    rmd_content <- c(rmd_content,
      "## Diversity Retention",
      ""
    )
    if (!is.na(qc_metrics$shannon_ens_retention_percent)) {
      rmd_content <- c(rmd_content, sprintf("- Shannon ENS (q=1): %.2f%% retained", qc_metrics$shannon_ens_retention_percent))
    }
    if (!is.na(qc_metrics$simpson_ens_retention_percent)) {
      rmd_content <- c(rmd_content, sprintf("- Simpson ENS (q=2): %.2f%% retained", qc_metrics$simpson_ens_retention_percent))
    }
    rmd_content <- c(rmd_content, "")
  }

  # Procrustes
  if (!is.na(qc_metrics$procrustes_m2)) {
    rmd_content <- c(rmd_content,
      "## Compositional Similarity",
      "",
      sprintf("- Procrustes M^2: %.6f *(lower = more similar)*", qc_metrics$procrustes_m2)
    )
    if (!is.na(qc_metrics$procrustes_correlation)) {
      rmd_content <- c(rmd_content, sprintf("- Procrustes Correlation: %.4f", qc_metrics$procrustes_correlation))
    }
    rmd_content <- c(rmd_content, "")
  }

  # Write Rmd and render
  writeLines(rmd_content, tmp_rmd)

  # Ensure output directory is absolute and exists
  output_dir_abs <- normalizePath(output_dir, mustWork = FALSE)
  if (!dir.exists(output_dir_abs)) {
    dir.create(output_dir_abs, recursive = TRUE)
  }

  pdf_path <- file.path(output_dir_abs, paste0(prefix, "_filtering_report.pdf"))

  tryCatch({
    # Save original working directory and set to output dir for rendering
    orig_wd <- getwd()
    on.exit(setwd(orig_wd), add = TRUE)
    setwd(output_dir_abs)
    rmarkdown::render(tmp_rmd, output_file = basename(pdf_path), quiet = !verbose)
    unlink(tmp_rmd)
    pdf_path
  }, error = function(e) {
    if (verbose) cat(sprintf("  Error generating PDF: %s\n", e$message))
    unlink(tmp_rmd)
    NULL
  })
}
