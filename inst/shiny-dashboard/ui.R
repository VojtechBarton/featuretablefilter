fluidPage(
  titlePanel("featuretablefilter - Interactive Filtering Dashboard"),

  # Include custom CSS
  includeCSS(system.file("shiny-dashboard/www/style.css", package = "featuretablefilter")),

  splitLayout(
    panelWidths = c(0.35, 0.65),

    # Sidebar: Controls and Parameters
    sidebarPanel(width = 12,

      h4("Filtering Parameters"),
      hr(),

      # Data Upload
      fileInput("data_file", "Upload Feature Table",
                accept = c(".tsv", ".csv", ".txt")),
      uiOutput("data_format_ui"),
      hr(),

      # Summary Statistics
      wellPanel(
        h5("Data Summary"),
        verbatimTextOutput("data_summary")
      ),
      hr(),

      # Coverage Filtering
      accordion(id = "coverage_accordion",

        accordionPanel("Coverage Filtering", icon = icon("filter"),

          selectInput("cov_method", "Method:",
                      choices = c("none", "absolute", "mad", "iqr", "good", "chao"),
                      selected = "mad"),

          uiOutput("cov_params_ui")
        ),

        accordionPanel("Singleton Ratio Filtering", icon = icon("exclamation-triangle"),

          selectInput("singleton_method", "Method:",
                      choices = c("none", "absolute"),
                      selected = "none"),

          uiOutput("singleton_params_ui")
        ),

        accordionPanel("Cross-Talk Filtering", icon = icon("exchange-alt"),

          selectInput("crosstalk_method", "Method:",
                      choices = c("none", "zero", "remove_feature", "flag"),
                      selected = "none"),

          uiOutput("crosstalk_params_ui")
        )
      ),
      hr(),

      # Abundance Filtering
      accordion(id = "abundance_accordion",

        accordionPanel("Abundance Filtering", icon = icon("arrow-down"),

          selectInput("abun_method", "Method:",
                      choices = c("none", "absolute", "relative", "relative_cutoff", "joint"),
                      selected = "joint"),

          uiOutput("abun_params_ui")
        )
      ),
      hr(),

      # Action Buttons
      actionButton("apply_filter", "Apply Filtering", class = "btn-primary btn-lg"),
      width = 4

    ),

    # Main Panel: Results and Visualizations
    mainPanel(width = 8,

      # Tabs for different views
      tabsetPanel(type = "tabs",

        # Tab 1: Overview
        tabPanel("Overview",
          fluidRow(
            column(6,
              wellPanel(
                h5("Original Table"),
                verbatimTextOutput("original_stats")
              )
            ),
            column(6,
              wellPanel(
                h5("Filtered Table"),
                verbatimTextOutput("filtered_stats")
              )
            )
          ),
          fluidRow(
            column(12,
              plotOutput("retention_plot", height = "300px")
            )
          )
        ),

        # Tab 2: Coverage
        tabPanel("Coverage Distribution",
          plotOutput("coverage_histogram", height = "400px")
        ),

        # Tab 3: Abundance
        tabPanel("Feature Abundance",
          plotOutput("abundance_histogram", height = "400px")
        ),

        # Tab 4: Sparsity
        tabPanel("Sparsity",
          plotOutput("sparsity_plot", height = "400px")
        ),

        # Tab 5: Scree Analysis
        tabPanel("Scree Analysis",
          plotOutput("scree_plot", height = "400px")
        ),

        # Tab 6: Results
        tabPanel("Results",
          fluidRow(
            column(12,
              h4("Filtered Feature Table Preview"),
              DT::dataTableOutput("filtered_table_preview")
            )
          ),
          hr(),
          fluidRow(
            column(4,
              downloadButton("download_table", "Download Table (TSV)")
            ),
            column(4,
              downloadButton("download_report", "Download Report")
            ),
            column(4,
              downloadButton("download_plots", "Download Plots (ZIP)")
            )
          )
        ),

        # Tab 7: R Code
        tabPanel("R Code",
          h4("Reproducible R Code"),
          p("Use this code to reproduce the filtering in your own script:"),
          verbatimTextOutput("r_code")
        )
      )
    )
  )
)
