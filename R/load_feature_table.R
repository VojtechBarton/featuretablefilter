#' Load a feature table from TSV or CSV file
#'
#' Automatically detects the file format and structure of microbiome feature tables.
#' Handles common formats where the first column contains feature IDs.
#'
#' @param file Path to the input file (TSV or CSV).
#' @param header Logical. Should the first row be treated as headers (sample names)?
#'               If NULL, will be auto-detected.
#' @param sep File separator. If NULL, will be auto-detected (.tsv uses tab, .csv uses comma).
#' @param feature_col Name or index of the column containing feature IDs. Default is 1 (first column).
#' @param feature_name_pattern Optional regex pattern to identify the feature ID column by name.
#'                             Common patterns: "OTU", "feature", "ASV", "taxon".
#'
#' @return A data.frame with the feature table. Feature IDs are kept as a regular column
#'         (not row names) to preserve the standard OTU/feature table format.
#'
#' @export
#'
#' @examples
#' # Load TSV file (auto-detects format)
#' # table <- load_feature_table("example_feature_table.tsv")
#'
#' # Load CSV file
#' # table <- load_feature_table("feature-table.csv")
#'
#' # Specify feature column by name
#' # table <- load_feature_table("table.tsv", feature_col = "OTU_ID")
load_feature_table <- function(file, header = NULL, sep = NULL,
                                feature_col = 1, feature_name_pattern = NULL) {
  # Check file exists
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  # Auto-detect separator based on file extension
  if (is.null(sep)) {
    ext <- tolower(tools::file_ext(file))
    sep <- if (ext == "tsv" || ext == "tab") "\t" else ","
  }

  # Read first line to auto-detect header
  if (is.null(header)) {
    first_line <- readLines(file, n = 1, warn = FALSE)
    first_fields <- strsplit(first_line, sep)[[1]]

    # Check if first field looks like a feature ID (not numeric)
    first_field <- first_fields[1]
    is_feature_id <- !grepl("^[0-9.e+-]+$", trimws(first_field))

    # Check other fields - if they look like sample names (not all numeric), it has header
    other_fields <- first_fields[-1]
    non_numeric_samples <- sum(!grepl("^[0-9.e+-]+$", trimws(other_fields)))

    header <- is_feature_id && non_numeric_samples > length(other_fields) / 2
  }

  # Read the table
  # Note: comment.char = "" is critical because feature IDs may start with # (e.g., #OTU ID)
  table_data <- read.table(
    file,
    header = header,
    sep = sep,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""
  )

  # Try to detect feature column by pattern if not specified by position
  if (is.character(feature_col) && length(feature_col) == 1) {
    if (feature_col %in% colnames(table_data)) {
      # Column name found
    } else if (!is.null(feature_name_pattern)) {
      matched <- grep(feature_name_pattern, colnames(table_data), value = TRUE, ignore.case = TRUE)
      if (length(matched) > 0) {
        feature_col <- matched[1]
      }
    } else {
      stop("Feature column '", feature_col, "' not found in table")
    }
  }

  # Ensure feature column is first
  if (is.numeric(feature_col) && feature_col != 1) {
    cols <- c(feature_col, setdiff(1:ncol(table_data), feature_col))
    table_data <- table_data[, cols]
  } else if (is.character(feature_col) && colnames(table_data)[1] != feature_col) {
    cols <- c(feature_col, setdiff(colnames(table_data), feature_col))
    table_data <- table_data[, cols]
  }

  return(table_data)
}
