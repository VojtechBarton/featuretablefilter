#' Convert absolute abundance table to relative abundance table
#'
#' Transforms a feature table from absolute counts (reads) to relative abundances
#' (proportions) by normalizing each sample to sum to 1 (or 100 if percent = TRUE).
#'
#' @param table A feature table (data.frame or matrix) with features as rows and samples as columns.
#'              First column should be feature IDs, remaining columns are sample counts.
#' @param percent Logical. If TRUE, return percentages (sum to 100). If FALSE (default),
#'                return proportions (sum to 1).
#' @param zero_replace Value to replace NA/NaN values that may arise from dividing 0/0.
#'                     Default is 0.
#'
#' @return A data.frame with the same structure as input, but with relative abundances
#'         instead of absolute counts.
#'
#' @export
#'
#' @examples
#' # Convert to proportions
#' # rel_table <- convert_to_relative(table)
#'
#' # Convert to percentages
#' # rel_table <- convert_to_relative(table, percent = TRUE)
convert_to_relative <- function(table, percent = FALSE, zero_replace = 0) {
  # Extract feature IDs and count columns
  feature_ids <- table[, 1, drop = FALSE]
  counts <- table[, -1, drop = FALSE]

  # Calculate sample totals
  sample_totals <- colSums(counts)

  # Convert to relative abundances
  # Use sweep for efficient column-wise division
  relative <- sweep(counts, 2, sample_totals, FUN = "/")

  # Replace NA/NaN (from 0/0 when sample has no reads)
  relative[is.na(relative)] <- zero_replace

  # Convert to percentages if requested
  if (percent) {
    relative <- relative * 100
  }

  # Combine feature IDs with relative abundances
  result <- cbind(feature_ids, relative)

  return(result)
}
