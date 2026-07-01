library(featuretablefilter)
library(DT)

server <- function(input, output, session) {

  # Reactive values
  rv <- reactiveValues(
    original_table = NULL,
    filtered_table = NULL,
    qc_metrics = NULL,
    filtering_params = NULL
  )

  # Data upload handling
  output$data_format_ui <- renderUI({
    req(input$data_file)
    selectInput("data_format", "File Format:",
                choices = c("auto", "tsv", "csv"),
                selected = "auto")
  })

  # Load data when file is uploaded
  data <- reactive({
    req(input$data_file)

    format <- input$data_format
    if (format == "auto" || format == "tsv") {
      table <- load_feature_table(input$data_file$datapath, sep = "\t")
    } else {
      table <- load_feature_table(input$data_file$datapath, sep = ",")
    }
    table
  })

  # Update data summary
  output$data_summary <- renderText({
    req(data())
    tbl <- data()
    paste0(
      "Features: ", nrow(tbl), "\n",
      "Samples: ", ncol(tbl) - 1, "\n",
      "Total reads: ", format(sum(tbl[, -1]), big.mark = ","), "\n",
      "Sparsity: ", sprintf("%.1f%%", mean(tbl[, -1] == 0) * 100)
    )
  })

  # Store original table when data loads
  observe({
    req(data())
    rv$original_table <- data()
    rv$filtered_table <- data()
  })

  # Dynamic UI for coverage parameters
  output$cov_params_ui <- renderUI({
    req(rv$original_table)

    method <- input$cov_method
    if (method == "none") {
      return(NULL)
    }

    tagList(
      if (method == "absolute" || method == "mad" || method == "iqr") {
        numericInput("cov_threshold",
                     if (method == "absolute") "Min reads:" else "Multiplier:",
                     value = if (method == "absolute") 1000 else 3,
                     min = 0)
      },
      if (method == "mad" || method == "iqr") {
        numericInput("cov_floor", "Floor (min cutoff):", value = 0, min = 0)
      },
      if (method == "good" || method == "chao") {
        sliderInput("cov_target", "Target coverage:",
                    min = 0.5, max = 0.99, value = 0.95, step = 0.01)
      },
      numericInput("cov_min_reads", "Additional min reads floor:",
                   value = 0, min = 0)
    )
  })

  # Dynamic UI for singleton parameters
  output$singleton_params_ui <- renderUI({
    req(rv$original_table)

    method <- input$singleton_method
    if (method == "none") return(NULL)

    tagList(
      numericInput("singleton_max_ratio", "Max ratio:", value = 0.1,
                   min = 0, max = 1, step = 0.01),
      radioButtons("singleton_count_type", "Count type:",
                   choices = c("both", "singleton", "doubleton"),
                   selected = "both")
    )
  })

  # Dynamic UI for cross-talk parameters
  output$crosstalk_params_ui <- renderUI({
    req(rv$original_table)

    method <- input$crosstalk_method
    if (method == "none") return(NULL)

    tagList(
      numericInput("crosstalk_threshold", "Relative threshold:",
                   value = 0.001, min = 0, max = 1, step = 0.0001),
      numericInput("crosstalk_min_abs", "Min absolute cutoff:",
                   value = NULL, min = 0)
    )
  })

  # Dynamic UI for abundance parameters
  output$abun_params_ui <- renderUI({
    req(rv$original_table)

    method <- input$abun_method
    if (method == "none") return(NULL)

    tagList(
      if (method == "absolute") {
        numericInput("abun_threshold", "Min reads:", value = 5, min = 0)
      },
      if (method == "relative") {
        numericInput("abun_threshold", "Relative abundance:",
                     value = 0.001, min = 0, max = 1, step = 0.0001)
      },
      if (method == "relative_cutoff") {
        tagList(
          numericInput("abun_threshold", "Relative threshold:",
                       value = 0.01, min = 0, max = 1, step = 0.001),
          numericInput("min_coverage_for_relative", "Min coverage sample:",
                       value = 1000, min = 0)
        )
      },
      if (method == "joint") {
        tagList(
          numericInput("abun_threshold", "Abundance threshold:",
                       value = 0.001, min = 0, max = 1, step = 0.0001),
          numericInput("abun_prevalence", "Prevalence threshold:",
                       value = 0.2, min = 0, max = 1, step = 0.01),
          radioButtons("abun_logic", "Logic:",
                       choices = c("OR", "AND"),
                       selected = "OR")
        )
      },
      numericInput("abun_min_samples", "Min samples:", value = 1, min = 1)
    )
  })

  # Apply filtering
  observeEvent(input$apply_filter, {
    req(rv$original_table)

    tryCatch({
      withProgress(message = "Applying filtering...", value = 0, {

        table <- rv$original_table

        # Coverage filtering
        if (input$cov_method != "none") {
          incProgress(0.2, detail = "Coverage filtering...")

          if (input$cov_method == "absolute") {
            table <- filter_by_coverage(table, min_reads = input$cov_threshold)
          } else if (input$cov_method == "mad") {
            est <- estimate_mad_cutoff(table,
                                       multiplier = input$cov_threshold,
                                       floor = input$cov_floor)
            table <- filter_by_coverage(table, min_reads = est$cutoff)
          } else if (input$cov_method == "iqr") {
            est <- estimate_iqr_cutoff(table,
                                       multiplier = input$cov_threshold,
                                       floor = input$cov_floor)
            table <- filter_by_coverage(table, min_reads = est$cutoff)
          } else if (input$cov_method %in% c("good", "chao")) {
            result <- filter_by_coverage_estimator(
              table,
              method = input$cov_method,
              target_coverage = input$cov_target,
              min_reads = input$cov_min_reads
            )
            table <- result$table
          }
        }

        # Singleton filtering
        if (input$singleton_method != "none") {
          incProgress(0.2, detail = "Singleton filtering...")
          table <- filter_by_singleton_ratio(
            table,
            max_singleton_ratio = input$singleton_max_ratio,
            count_type = input$singleton_count_type
          )
        }

        # Cross-talk filtering
        if (input$crosstalk_method != "none") {
          incProgress(0.2, detail = "Cross-talk filtering...")
          table <- filter_cross_talk(
            table,
            max_rel_threshold = input$crosstalk_threshold,
            min_abs_cutoff = input$crosstalk_min_abs,
            mode = input$crosstalk_method
          )
        }

        # Abundance filtering
        if (input$abun_method != "none") {
          incProgress(0.3, detail = "Abundance filtering...")

          if (input$abun_method == "absolute") {
            table <- filter_features_by_abundance(
              table,
              threshold = input$abun_threshold,
              mode = "absolute",
              min_samples = input$abun_min_samples
            )
          } else if (input$abun_method == "relative") {
            table <- filter_features_by_abundance(
              table,
              threshold = input$abun_threshold,
              mode = "relative",
              min_samples = input$abun_min_samples
            )
          } else if (input$abun_method == "relative_cutoff") {
            result <- filter_by_relative_cutoff(
              table,
              min_coverage = input$min_coverage_for_relative,
              relative_threshold = input$abun_threshold
            )
            table <- result$table
          } else if (input$abun_method == "joint") {
            result <- filter_features_joint(
              table,
              abundance_threshold = input$abun_threshold,
              prevalence_threshold = input$abun_prevalence,
              mode = "relative",
              logic = input$abun_logic
            )
            table <- result$table
          }
        }

        incProgress(0.1, detail = "Computing QC metrics...")
        qc <- compute_filtering_qc(rv$original_table, table)

        rv$filtered_table <- table
        rv$qc_metrics <- qc

        # Store params for code generation
        rv$filtering_params <- list(
          cov_method = input$cov_method,
          cov_threshold = input$cov_threshold,
          abun_method = input$abun_method,
          abun_threshold = input$abun_threshold
        )
      })
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })

  # Original stats
  output$original_stats <- renderText({
    req(rv$original_table)
    tbl <- rv$original_table
    paste0(
      "Features: ", nrow(tbl), "\n",
      "Samples: ", ncol(tbl) - 1, "\n",
      "Total reads: ", format(sum(tbl[, -1]), big.mark = ",")
    )
  })

  # Filtered stats
  output$filtered_stats <- renderText({
    req(rv$filtered_table)
    tbl <- rv$filtered_table
    paste0(
      "Features: ", nrow(tbl), "\n",
      "Samples: ", ncol(tbl) - 1, "\n",
      "Total reads: ", format(sum(tbl[, -1]), big.mark = ",")
    )
  })

  # Retention plot
  output$retention_plot <- renderPlot({
    req(rv$qc_metrics)
    qc <- rv$qc_metrics

    df <- data.frame(
      Metric = c("Features", "Samples", "Reads"),
      Retention = c(
        qc$feature_retention_percent,
        qc$sample_retention_percent,
        qc$read_retention_percent
      )
    )

    ggplot(df, aes(x = Metric, y = Retention, fill = Metric)) +
      geom_col() +
      geom_text(aes(label = sprintf("%.1f%%", Retention)), vjust = -0.5) +
      scale_y_continuous(limits = c(0, 100)) +
      theme_minimal() +
      labs(title = "Retention Rates", y = "Percentage (%)") +
      theme(legend.position = "none")
  })

  # Coverage histogram
  output$coverage_histogram <- renderPlot({
    req(rv$original_table)
    plot_coverage_histogram(rv$original_table)
  })

  # Abundance histogram
  output$abundance_histogram <- renderPlot({
    req(rv$original_table)
    tbl <- rv$original_table
    abundances <- rowSums(tbl[, -1])

    ggplot(data.frame(abundances), aes(x = abundances)) +
      geom_histogram(binwidth = 10, fill = "steelblue", color = "white") +
      scale_x_log10() +
      theme_minimal() +
      labs(title = "Feature Abundance Distribution",
           x = "Total Reads (log10)", y = "Frequency")
  })

  # Sparsity plot
  output$sparsity_plot <- renderPlot({
    req(rv$original_table, rv$filtered_table)

    orig_sparsity <- apply(rv$original_table[, -1], 2, function(x) mean(x == 0))
    filt_sparsity <- apply(rv$filtered_table[, -1], 2, function(x) mean(x == 0))

    df <- data.frame(
      Sparsity = c(orig_sparsity, filt_sparsity),
      Group = rep(c("Original", "Filtered"),
                  each = length(orig_sparsity))
    )

    ggplot(df, aes(x = Sparsity, fill = Group)) +
      geom_density(alpha = 0.5) +
      theme_minimal() +
      labs(title = "Sparsity Distribution",
           x = "Proportion Zeros", y = "Density")
  })

  # Scree plot
  output$scree_plot <- renderPlot({
    req(rv$filtered_table)

    scree <- compute_scree(rv$filtered_table, type = "mad_multiplier", n_steps = 20)
    plot_scree(scree)
  })

  # Filtered table preview
  output$filtered_table_preview <- DT::renderDataTable({
    req(rv$filtered_table)
    DT::datatable(head(rv$filtered_table, 100), options = list(pageLength = 10))
  })

  # Generate R code
  output$r_code <- renderText({
    req(rv$filtering_params)
    p <- rv$filtering_params

    code <- 'library(featuretablefilter)

# Load data
table <- load_feature_table("your_data.tsv")

# Apply filtering
result <- run_filtering_pipeline(
  input = table,'

    if (p$cov_method != "none") {
      code <- paste0(code, "\n  cov_filter_method = \"", p$cov_method, "\",")
    }
    if (p$abun_method != "none") {
      code <- paste0(code, "\n  abun_filter_method = \"", p$abun_method, "\",")
    }

    code <- paste0(code, '
  generate_plots = TRUE,
  generate_report = TRUE
)

# Access results
filtered <- result$filtered_table
')
    code
  })

  # Downloads
  output$download_table <- downloadHandler(
    filename = function() {"filtered_table.tsv"},
    content = function(file) {
      req(rv$filtered_table)
      write.table(rv$filtered_table, file, sep = "\t",
                  row.names = FALSE, quote = FALSE)
    }
  )

  output$download_report <- downloadHandler(
    filename = function() {"filtering_report.txt"},
    content = function(file) {
      req(rv$qc_metrics)
      lines <- c(
        "Filtering Report",
        "================",
        "",
        "QC Metrics:",
        paste("  Feature retention:", sprintf("%.1f%%", rv$qc_metrics$feature_retention_percent)),
        paste("  Sample retention:", sprintf("%.1f%%", rv$qc_metrics$sample_retention_percent)),
        paste("  Read retention:", sprintf("%.1f%%", rv$qc_metrics$read_retention_percent))
      )
      writeLines(lines, file)
    }
  )

  output$download_plots <- downloadHandler(
    filename = function() {"plots.zip"},
    content = function(file) {
      # Placeholder - would need to generate and zip plots
      showNotification("Plots download coming soon!", type = "message")
    }
  )
}
