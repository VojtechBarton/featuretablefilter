#' Run the featuretablefilter Interactive Dashboard
#'
#' Launches a Shiny web application for interactive feature table filtering.
#' The dashboard allows you to explore different filtering methods and parameters
#' with real-time preview of results and visualizations.
#'
#' @param host Host address to run the app on. Default is "127.0.0.1" (localhost).
#' @param port Port number to run the app on. Default is 0 (auto-select available port).
#' @param launch.browser Logical. Whether to automatically launch the browser.
#'                         Default is TRUE.
#' @param max_upload_size Maximum file upload size in bytes. Default is 100 MB
#'                        (100 * 1024 * 1024 = 104857600 bytes). Set higher for
#'                        large feature tables.
#' @param ... Additional arguments passed to [shiny::runApp()].
#'
#' @return NULL (invisibly)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Launch the dashboard in default browser (100 MB max upload)
#' runDashboard()
#'
#' # Launch on specific port
#' runDashboard(port = 3838)
#'
#' # Launch without auto-opening browser
#' runDashboard(launch.browser = FALSE)
#'
#' # Launch with larger upload limit (500 MB)
#' runDashboard(max_upload_size = 500 * 1024 * 1024)
#' }
runDashboard <- function(host = "127.0.0.1", port = 0,
                          launch.browser = TRUE, max_upload_size = 100 * 1024 * 1024,
                          ...) {
  # Check for shiny
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("shiny package is required. Install with: install.packages('shiny')")
  }

  # Set maximum upload size (default: 100 MB)
  old_option <- options(shiny.maxRequestSize = max_upload_size)
  on.exit(options(old_option), add = TRUE)

  # Get the app directory
  app_dir <- system.file("shiny-dashboard", package = "featuretablefilter")

  if (app_dir == "") {
    stop(
      "Could not find shiny-dashboard directory. ",
      "The package may not have been installed correctly."
    )
  }

  # Run the Shiny app
  shiny::runApp(
    appDir = app_dir,
    host = host,
    port = port,
    launch.browser = launch.browser,
    ...
  )
}
