#' Calculate sparsity of a feature table
#'
#' Computes the proportion of zero cells in the table, which is a common
#' quality control metric for microbiome and other sparse count data.
#'
#' @param table A feature table (data.frame or matrix) with features as rows and samples as columns.
#'              First column should be feature IDs, remaining columns are sample counts/abundances.
#' @param by_sample Logical. If TRUE, returns sparsity per sample. If FALSE (default),
#'                  returns overall table sparsity.
#'
#' @return If by_sample = FALSE: a numeric value between 0 and 1 representing
#'         the proportion of zero cells in the table.
#'         If by_sample = TRUE: a named numeric vector with sparsity for each sample.
#'
#' @export
#'
#' @examples
#' # Overall table sparsity
#' # sparsity <- calculate_sparsity(table)
#'
#' # Sparsity per sample
#' # sample_sparsity <- calculate_sparsity(table, by_sample = TRUE)
calculate_sparsity <- function(table, by_sample = FALSE) {
  # Extract abundance columns (skip first column which is feature IDs)
  abundances <- table[, -1, drop = FALSE]

  if (by_sample) {
    # Calculate sparsity per sample
    n_zeros <- colSums(abundances == 0)
    n_total <- nrow(abundances)
    sparsity <- n_zeros / n_total
    names(sparsity) <- colnames(abundances)
    return(sparsity)
  } else {
    # Calculate overall sparsity
    n_zeros <- sum(abundances == 0)
    n_total <- length(abundances)
    return(n_zeros / n_total)
  }
}
