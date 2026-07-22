#' Calculate coefficient of variation (CV) for features
#'
#' Computes the coefficient of variation (standard deviation / mean) for each feature
#' across samples. This is a useful metric for identifying highly variable features.
#'
#' @param table A feature table (data.frame or matrix) with features as rows and samples as columns.
#'              First column should be feature IDs, remaining columns are sample counts/abundances.
#' @param plot Logical. If TRUE, creates a histogram of CV distribution. Default is FALSE.
#' @param plot_type Type of plot: "histogram" (default) or "boxplot".
#' @param log_transform Logical. If TRUE, applies log10 transformation to CV values for plotting.
#'                      Useful when CV values span several orders of magnitude. Default is FALSE.
#' @param main Plot title (only used if plot = TRUE).
#'
#' @return A data.frame with feature IDs and their corresponding CV values.
#'         Columns: feature_id, mean, sd, cv
#'
#' @export
#'
#' @examples
#' data(example_feature_table)
#' result <- calculate_feature_cv(example_feature_table)
#' head(result)
calculate_feature_cv <- function(table, plot = FALSE, plot_type = c("histogram", "boxplot"),
                                  log_transform = FALSE, main = "Feature Coefficient of Variation") {
  # Extract feature IDs and abundance columns
  feature_ids <- table[, 1]
  abundances <- as.matrix(table[, -1, drop = FALSE])

  # Calculate mean and SD for each feature
  feature_means <- rowMeans(abundances)
  feature_sds <- apply(abundances, 1, sd)

  # Calculate CV (avoid division by zero)
  cv <- feature_sds / feature_means
  cv[is.infinite(cv)] <- NA
  cv[is.na(cv)] <- 0

  # Create result data.frame
  result <- data.frame(
    feature_id = feature_ids,
    mean = feature_means,
    sd = feature_sds,
    cv = cv,
    stringsAsFactors = FALSE
  )

  # Optional plotting
  if (plot) {
    plot_type <- match.arg(plot_type)

    cv_values <- result$cv
    if (log_transform) {
      # Add small constant to avoid log(0)
      cv_values <- log10(cv_values + 1e-10)
      xlab <- "log10(CV + 1e-10)"
      main <- paste(main, "(log10 scale)")
    } else {
      xlab <- "Coefficient of Variation"
    }

    if (plot_type == "histogram") {
      hist(cv_values,
           main = main,
           xlab = xlab,
           col = "steelblue",
           border = "white",
           breaks = 30)
    } else {
      boxplot(cv_values,
              main = main,
              ylab = "Coefficient of Variation",
              col = "steelblue",
              border = "white")
    }
  }

  return(result)
}
