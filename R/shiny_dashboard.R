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
#' @param ... Additional arguments passed to [shiny::runApp()].
#'
#' @return NULL (invisibly)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Launch the dashboard in default browser
#' runDashboard()
#'
#' # Launch on specific port
#' runDashboard(port = 3838)
#'
#' # Launch without auto-opening browser
#' runDashboard(launch.browser = FALSE)
#' }
runDashboard <- function(host = "127.0.0.1", port = 0,
                          launch.browser = TRUE, ...) {
  # Check for shiny
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("shiny package is required. Install with: install.packages('shiny')")
  }

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
