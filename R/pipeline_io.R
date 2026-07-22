# =============================================================================
# Pipeline Input/Output Helper Functions (Internal)
# =============================================================================
# These functions handle loading input data and saving pipeline outputs.
# All functions are internal (prefixed with '.') and not exported.
# =============================================================================

#' Load pipeline input data from file or object#'#' Handles multiple input types: file paths, phyloseq objects, TreeSummarizedExperiment,#' data.frame, or matrix. Returns a list with the table and input metadata.#'#' @param input Input source (file path string, phyloseq, TSE, data.frame, or matrix)#' @param verbose Logical. Print progress messages?#'#' @return A list with components:#' \describe{#'   \item{table}{Feature table as data.frame}#'   \item{input_class}{Original input class ("data.frame", "phyloseq", "TreeSummarizedExperiment")}#'   \item{input_file_path}{File path if input was a file, NULL otherwise}#'   \item{original_object}{Original input object (for later conversion back)}#' }
#' @noRd
.load_pipeline_input <- function(input, verbose = TRUE) {
  input_class <- NULL
  input_file_path <- NULL
  original_object <- NULL

  if (is.character(input) && length(input) == 1 && file.exists(input)) {
    # Input is a file path
    input_class <- "data.frame"
    if (verbose) cat(sprintf("Loading feature table from file: %s\n", input))
    table_data <- load_feature_table(input)
    input_file_path <- input

  } else if (inherits(input, "phyloseq")) {
    # Input is a phyloseq object
    input_class <- "phyloseq"
    if (verbose) cat("Converting phyloseq object to feature table...\n")
    table_data <- from_phyloseq(input, include_taxa = FALSE)
    original_object <- input

  } else if (inherits(input, c("TreeSummarizedExperiment", "SingleCellExperiment", "SummarizedExperiment"))) {
    # Input is a TreeSummarizedExperiment or related object
    input_class <- "TreeSummarizedExperiment"
    if (verbose) cat("Converting TreeSummarizedExperiment object to feature table...\n")
    table_data <- from_TSE(input, add_row_data = FALSE)
    original_object <- input

  } else if (is.data.frame(input) || is.matrix(input)) {
    # Input is already a data.frame or matrix
    input_class <- "data.frame"
    if (verbose) cat("Using provided data.frame/matrix as feature table...\n")
    table_data <- as.data.frame(input)
    if (is.matrix(input)) {
      feature_ids <- rownames(table_data)
      if (is.null(feature_ids)) {
        feature_ids <- paste0("Feature_", seq_len(nrow(table_data)))
      }
      table_data <- data.frame(feature_id = feature_ids, table_data, stringsAsFactors = FALSE)
    }

  } else {
    stop("Input must be a file path (character), data.frame, matrix, phyloseq, or TreeSummarizedExperiment object")
  }

  list(
    table = table_data,
    input_class = input_class,
    input_file_path = input_file_path,
    original_object = original_object
  )
}

#' Create output directory if it doesn't exist#'#' Only creates directory for file-based inputs (not for object inputs).#'#' @param output_dir Directory path to create#' @param input_file_path File path from input (NULL if input was an object)#' @param verbose Logical. Print progress messages?#'#' @return The output directory path (invariant)
#' @noRd
.create_output_directory <- function(output_dir, input_file_path, verbose = TRUE) {
  # Only create directory for file inputs
  if (!is.null(input_file_path)) {
    if (!dir.exists(output_dir)) {
      if (verbose) cat(sprintf("Creating output directory: %s\n", output_dir))
      dir.create(output_dir, recursive = TRUE)
    }
  }
  output_dir
}

#' Save filtered table to file#'#' Writes the filtered feature table to a TSV file.#'#' @param table Filtered feature table (data.frame)#' @param output_dir Output directory path#' @param prefix Prefix for output filename#' @param verbose Logical. Print progress messages?#'#' @return Path to the saved file
#' @noRd
.save_filtered_table <- function(table, output_dir, prefix, verbose = TRUE) {
  output_path <- file.path(output_dir, paste0(prefix, "_filtered_table.tsv"))
  write.table(table, output_path, sep = "\t", row.names = FALSE, quote = FALSE)
  if (verbose) cat(sprintf("Filtered table saved to: %s\n", output_path))
  output_path
}

#' Convert filtered table back to original input class#'#' Converts the filtered data.frame back to phyloseq or TreeSummarizedExperiment#' if the original input was one of those classes.#'#' @param filtered_table Filtered feature table (data.frame)#' @param original_class Original input class#' @param original_object Original input object (phyloseq or TSE)#'#' @return A list with original_table and filtered_table in original class format
#' @noRd
.convert_output_back_to_original_class <- function(filtered_table, original_class, original_object) {
  original_table_out <- filtered_table
  filtered_table_out <- filtered_table

  if (original_class == "phyloseq") {
    original_table_out <- original_object
    filtered_table_out <- to_phyloseq(
      filtered_table,
      tax_table = phyloseq::tax_table(original_object),
      phy_tree = phyloseq::phy_tree(original_object),
      sample_data = phyloseq::sample_data(original_object)
    )

  } else if (original_class == "TreeSummarizedExperiment") {
    original_table_out <- original_object
    filtered_table_out <- to_TSE(
      filtered_table,
      rowData = SummarizedExperiment::rowData(original_object),
      colData = SummarizedExperiment::colData(original_object),
      reducedDims = SummarizedExperiment::reducedDims(original_object),
      rowTree = TreeSummarizedExperiment::rowTree(original_object),
      rowLinks = TreeSummarizedExperiment::rowLinks(original_object)
    )
  }

  list(
    original_table = original_table_out,
    filtered_table = filtered_table_out
  )
}

#' Generate filtering summary step entry#'#' Creates a standardized summary entry for a filtering step.#'#' @param step_name Name of the filtering step#' @param method Method used for filtering#' @param params Named list of parameters used#' @param before List with samples_before, features_before, reads_before#' @param after List with samples_after, features_after, reads_after#'#' @return A named list with step summary information
#' @noRd
.format_step_summary <- function(step_name, method, params, before, after) {
  list(
    step = step_name,
    method = method,
    params = params,
    samples_before = before$samples,
    samples_after = after$samples,
    features_before = before$features,
    features_after = after$features,
    reads_before = before$reads,
    reads_after = after$reads,
    samples_removed = before$samples - after$samples,
    features_removed = before$features - after$features,
    reads_removed = before$reads - after$reads
  )
}
