#' Convert from phyloseq object to feature table data.frame
#'
#' Extracts the OTU table from a phyloseq object and converts it to the standard
#' data.frame format used by this package (feature IDs in first column, samples as columns).
#'
#' @param phylo_obj A phyloseq object containing at least an OTU table.
#' @param transpose Logical. If TRUE (default), transposes so features are rows and samples are columns.
#'                  phyloseq stores OTU tables with features as rows by default, but some operations
#'                  may change this.
#' @param include_taxa Logical. If TRUE and tax_table is present in phyloseq, appends taxonomy
#'                     information as additional columns after the feature ID.
#'
#' @return A data.frame with feature IDs in the first column, followed by sample counts.
#'         If include_taxa=TRUE, taxonomy columns are appended.
#'
#' @export
#'
#' @examples
#' # if (requireNamespace("phyloseq", quietly = TRUE)) {
#' #   ps <- phyloseq::phyloseq(
#' #     phyloseq::otu_table(matrix(rpois(100, 5), 10, 10), taxa_are_rows = TRUE),
#' #     phyloseq::tax_table(matrix(sample(c("A","B","C"), 100, replace=TRUE), 10, 10))
#' #   )
#' #   df <- from_phyloseq(ps)
#' # }
from_phyloseq <- function(phylo_obj, transpose = TRUE, include_taxa = FALSE) {
  # Check for phyloseq
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("phyloseq package is required. Install with: install.packages('phyloseq')")
  }

  # Validate input
  if (!phyloseq::is.phyloseq(phylo_obj)) {
    stop("Input must be a phyloseq object")
  }

  # Extract OTU table
  otu <- phyloseq::otu_table(phylo_obj)
  if (is.null(otu)) {
    stop("phyloseq object must contain an OTU table")
  }

  # Convert to matrix if needed
  otu_mat <- as.matrix(otu)

  # Ensure features are rows, samples are columns
  if (phyloseq::taxa_are_rows(otu)) {
    # Already in correct orientation
    features <- rownames(otu_mat)
    samples <- colnames(otu_mat)
  } else {
    # Need to transpose
    otu_mat <- t(otu_mat)
    features <- rownames(otu_mat)
    samples <- colnames(otu_mat)
  }

  # Create feature ID column
  if (is.null(features)) {
    features <- paste0("Feature_", seq_len(nrow(otu_mat)))
  }

  result <- data.frame(
    feature_id = features,
    otu_mat,
    stringsAsFactors = FALSE
  )

  # Add taxonomy if requested and available
  if (include_taxa) {
    tax <- phyloseq::tax_table(phylo_obj)
    if (!is.null(tax)) {
      tax_mat <- as.matrix(tax)

      # Align taxonomy with OTU table
      if (phyloseq::taxa_are_rows(tax)) {
        tax_aligned <- tax_mat[features, , drop = FALSE]
      } else {
        tax_aligned <- tax_mat[, features, drop = FALSE]
        tax_aligned <- t(tax_aligned)
      }

      # Handle missing taxonomy entries
      colnames(tax_aligned) <- paste0("tax_", colnames(tax_aligned))
      tax_aligned[is.na(tax_aligned)] <- ""

      result <- cbind(result, tax_aligned)
    }
  }

  # Preserve sample names as column names
  if (!is.null(samples)) {
    colnames(result)[-1] <- samples
  }

  return(result)
}


#' Convert feature table data.frame to phyloseq object
#'
#' Creates a phyloseq object from a standard feature table data.frame.
#' Optionally incorporates taxonomy, phylogenetic tree, and sample data
#' if provided.
#'
#' @param table A data.frame with feature IDs in the first column and samples as columns.
#' @param tax_table Optional. Taxonomy table (data.frame or matrix) with feature IDs as rows
#'                  and taxonomic ranks as columns (Kingdom, Phylum, Class, Order, Family, Genus, Species).
#' @param phy_tree Optional. A phylogenetic tree (phylo object) with tip labels matching feature IDs.
#' @param sample_data Optional. Sample metadata (data.frame) with row names matching sample names.
#' @param feature_col Name or index of the column containing feature IDs. Default is "feature_id" or 1.
#'
#' @return A phyloseq object containing OTU table and optionally taxonomy, tree, and sample data.
#'
#' @export
#'
#' @examples
#' # Create simple phyloseq from feature table
#' # table <- data.frame(feature_id = paste0("ASV_", 1:10), matrix(rpois(100, 5), 10, 10))
#' # ps <- to_phyloseq(table)
#'
#' # With taxonomy
#' # tax <- data.frame(
#' #   feature_id = paste0("ASV_", 1:10),
#' #   Kingdom = rep("Bacteria", 10),
#' #   Phylum = sample(c("Firmicutes", "Bacteroidetes"), 10, replace = TRUE)
#' # )
#' # ps <- to_phyloseq(table, tax_table = tax)
to_phyloseq <- function(table, tax_table = NULL, phy_tree = NULL,
                         sample_data = NULL, feature_col = NULL) {
  # Check for phyloseq
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("phyloseq package is required. Install with: install.packages('phyloseq')")
  }

  # Validate input
  if (!is.data.frame(table) && !is.matrix(table)) {
    stop("table must be a data.frame or matrix")
  }

  table <- as.data.frame(table)

  # Detect feature column
  if (is.null(feature_col)) {
    if ("feature_id" %in% colnames(table)) {
      feature_col <- "feature_id"
    } else if (ncol(table) > 1 && !grepl("^[0-9.e+-]+$", colnames(table)[1])) {
      feature_col <- 1
    } else {
      stop("Could not detect feature ID column. Please specify feature_col.")
    }
  }

  # Extract feature IDs and OTU table
  if (is.numeric(feature_col)) {
    feature_ids <- table[, feature_col]
    otu_mat <- as.matrix(table[, -feature_col, drop = FALSE])
  } else {
    feature_ids <- table[[feature_col]]
    otu_mat <- as.matrix(table[, -which(colnames(table) == feature_col), drop = FALSE])
  }

  rownames(otu_mat) <- feature_ids

  # Create OTU table (with taxa_are_rows = TRUE for phyloseq standard)
  otu <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)

  # Build phyloseq components
  components <- list(otu_table = otu)

  # Add taxonomy if provided
  if (!is.null(tax_table)) {
    if (!requireNamespace("phyloseq", quietly = TRUE)) {
      warning("phyloseq required for taxonomy. Ignoring tax_table.")
    } else {
      tax_df <- as.data.frame(tax_table)

      # Check if tax_table has feature ID column
      if ("feature_id" %in% colnames(tax_df)) {
        tax_df <- tax_df[-which(colnames(tax_df) == "feature_id")]
      }

      # Align with OTU table
      tax_df <- tax_df[feature_ids, , drop = FALSE]
      rownames(tax_df) <- feature_ids

      # Convert to tax_table format
      tax_mat <- as.matrix(tax_df)
      tax <- phyloseq::tax_table(tax_mat)
      components$tax_table <- tax
    }
  }

  # Add phylogenetic tree if provided
  if (!is.null(phy_tree)) {
    if (!inherits(phy_tree, "phylo")) {
      warning("phy_tree must be a phylo object. Ignoring.")
    } else {
      components$phy_tree <- phy_tree
    }
  }

  # Add sample data if provided
  if (!is.null(sample_data)) {
    sd_df <- as.data.frame(sample_data)

    # Use row names if no explicit sample ID column
    if (is.null(rownames(sd_df))) {
      sample_names_in_otu <- colnames(otu_mat)
      if (!is.null(sample_names_in_otu)) {
        rownames(sd_df) <- sample_names_in_otu
      }
    }

    # Align with OTU samples
    sample_names_in_otu <- colnames(otu_mat)
    sd_aligned <- sd_df[sample_names_in_otu, , drop = FALSE]

    components$sam_data <- phyloseq::sample_data(sd_aligned)
  }

  # Create phyloseq object
  do.call(phyloseq::phyloseq, components)
}


#' Convert from TreeSummarizedExperiment to feature table data.frame
#'
#' Extracts the assay data from a TreeSummarizedExperiment object and converts it to
#' the standard data.frame format used by this package.
#'
#' @param tse A TreeSummarizedExperiment object.
#' @param assay_name Name of the assay to extract. Default is "counts" or the first available assay.
#' @param add_row_data Logical. If TRUE and rowData contains taxonomy, appends it as columns.
#' @param add_col_data Logical. If TRUE and colData contains sample metadata, returns a list
#'                     with both the feature table and sample data.
#'
#' @return A data.frame with feature IDs in the first column, followed by sample counts.
#'         If add_col_data=TRUE, returns a list with $table and $sample_data.
#'
#' @export
#'
#' @examples
#' # if (requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
#' #   # Assuming tse is a TreeSummarizedExperiment object
#' #   df <- from_TSE(tse)
#' # }
from_TSE <- function(tse, assay_name = NULL, add_row_data = FALSE, add_col_data = FALSE) {
  # Check for TreeSummarizedExperiment
  if (!requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
    stop("TreeSummarizedExperiment package is required. Install with: BiocManager::install('TreeSummarizedExperiment')")
  }

  # Also need SummarizedExperiment
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package is required. Install with: BiocManager::install('SummarizedExperiment')")
  }

  # Validate input
  if (!inherits(tse, "TreeSummarizedExperiment") &&
      !inherits(tse, "SingleCellExperiment") &&
      !inherits(tse, "SummarizedExperiment")) {
    stop("Input must be a TreeSummarizedExperiment, SingleCellExperiment, or SummarizedExperiment object")
  }

  # Select assay
  if (is.null(assay_name)) {
    assay_names <- names(SummarizedExperiment::assays(tse))
    if ("counts" %in% assay_names) {
      assay_name <- "counts"
    } else {
      assay_name <- assay_names[1]
    }
  }

  # Extract assay
  assay_data <- SummarizedExperiment::assay(tse, assay_name)

  # Ensure matrix format with features as rows
  if (!is.matrix(assay_data)) {
    assay_data <- as.matrix(assay_data)
  }

  # Get feature IDs from rownames or rowData
  feature_ids <- rownames(assay_data)
  if (is.null(feature_ids)) {
    rd <- SummarizedExperiment::rowData(tse)
    if (!is.null(rd) && ncol(rd) > 0) {
      feature_ids <- rownames(rd)
    }
  }
  if (is.null(feature_ids)) {
    feature_ids <- paste0("Feature_", seq_len(nrow(assay_data)))
  }

  # Create result data.frame
  result <- data.frame(
    feature_id = feature_ids,
    assay_data,
    stringsAsFactors = FALSE
  )

  # Add rowData (taxonomy) if requested
  if (add_row_data) {
    rd <- SummarizedExperiment::rowData(tse)
    if (!is.null(rd) && nrow(rd) > 0) {
      rd_df <- as.data.frame(rd)
      rd_df[is.na(rd_df)] <- ""

      # Align with assay rows
      rd_aligned <- rd_df[feature_ids, , drop = FALSE]
      colnames(rd_aligned) <- paste0("rowData_", colnames(rd_aligned))

      result <- cbind(result, rd_aligned)
    }
  }

  # Return sample data separately if requested
  if (add_col_data) {
    cd <- SummarizedExperiment::colData(tse)
    if (!is.null(cd) && nrow(cd) > 0) {
      sample_data <- as.data.frame(cd)
      return(list(
        table = result,
        sample_data = sample_data
      ))
    }
  }

  return(result)
}


#' Convert feature table data.frame to TreeSummarizedExperiment
#'
#' Creates a TreeSummarizedExperiment object from a standard feature table data.frame.
#' Optionally incorporates rowData (taxonomy), colData (sample metadata), and reduced dimensions.
#'
#' @param table A data.frame with feature IDs in the first column and samples as columns.
#' @param rowData Optional. Feature metadata (taxonomy) as data.frame with row names matching feature IDs.
#' @param colData Optional. Sample metadata as data.frame with row names matching sample names.
#' @param reducedDims Optional. Named list of matrices for reduced dimension representations.
#' @param rowTree Optional. Phylogenetic tree (phylo object) with tip labels matching feature IDs.
#' @param rowLinks Optional. Data frame linking features to tree nodes.
#' @param feature_col Name or index of the column containing feature IDs. Default is "feature_id" or 1.
#' @param assay_name Name for the assay. Default is "counts".
#'
#' @return A TreeSummarizedExperiment object with the assay data and optional metadata.
#'
#' @export
#'
#' @examples
#' # Create simple TSE from feature table
#' # table <- data.frame(feature_id = paste0("ASV_", 1:10), matrix(rpois(100, 5), 10, 10))
#' # tse <- to_TSE(table)
#'
#' # With rowData (taxonomy) and colData (sample metadata)
#' # rowData <- data.frame(
#' #   Kingdom = rep("Bacteria", 10),
#' #   Phylum = sample(c("Firmicutes", "Bacteroidetes"), 10, replace = TRUE)
#' # )
#' # rownames(rowData) <- paste0("ASV_", 1:10)
#' #
#' # colData <- data.frame(
#' #   Condition = rep(c("Control", "Treatment"), 5)
#' # )
#' # rownames(colData) <- colnames(table)[-1]
#' #
#' # tse <- to_TSE(table, rowData = rowData, colData = colData)
to_TSE <- function(table, rowData = NULL, colData = NULL, reducedDims = NULL,
                    rowTree = NULL, rowLinks = NULL,
                    feature_col = NULL, assay_name = "counts") {
  # Check for required packages
  if (!requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
    stop("TreeSummarizedExperiment package is required. Install with: BiocManager::install('TreeSummarizedExperiment')")
  }

  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package is required. Install with: BiocManager::install('SummarizedExperiment')")
  }

  # Validate input
  if (!is.data.frame(table) && !is.matrix(table)) {
    stop("table must be a data.frame or matrix")
  }

  table <- as.data.frame(table)

  # Detect feature column
  if (is.null(feature_col)) {
    if ("feature_id" %in% colnames(table)) {
      feature_col <- "feature_id"
    } else if (ncol(table) > 1 && !grepl("^[0-9.e+-]+$", colnames(table)[1])) {
      feature_col <- 1
    } else {
      stop("Could not detect feature ID column. Please specify feature_col.")
    }
  }

  # Extract feature IDs and assay data
  if (is.numeric(feature_col)) {
    feature_ids <- table[, feature_col]
    assay_data <- as.matrix(table[, -feature_col, drop = FALSE])
    sample_cols <- colnames(table)[-feature_col]
  } else {
    feature_ids <- table[[feature_col]]
    assay_idx <- which(colnames(table) == feature_col)
    assay_data <- as.matrix(table[, -assay_idx, drop = FALSE])
    sample_cols <- colnames(table)[-assay_idx]
  }

  rownames(assay_data) <- feature_ids
  colnames(assay_data) <- sample_cols

  # Build assays list
  assays_list <- list()
  assays_list[[assay_name]] <- assay_data

  # Prepare rowData
  if (!is.null(rowData)) {
    rd_df <- as.data.frame(rowData)

    # Remove feature_id column if present
    if ("feature_id" %in% colnames(rd_df)) {
      rd_df <- rd_df[, -which(colnames(rd_df) == "feature_id"), drop = FALSE]
    }

    # Align with assay rows
    rd_aligned <- rd_df[feature_ids, , drop = FALSE]
    rownames(rd_aligned) <- feature_ids
    rowData <- S4Vectors::DataFrame(rd_aligned)
  }

  # Prepare colData
  if (!is.null(colData)) {
    cd_df <- as.data.frame(colData)

    # Use row names if not set
    if (is.null(rownames(cd_df))) {
      sample_names <- colnames(assay_data)
      if (!is.null(sample_names)) {
        rownames(cd_df) <- sample_names
      }
    }

    # Align with assay columns
    sample_names <- colnames(assay_data)
    cd_aligned <- cd_df[sample_names, , drop = FALSE]
    colData <- S4Vectors::DataFrame(cd_aligned)
  }

  # Handle reducedDims
  if (!is.null(reducedDims)) {
    if (!is.list(reducedDims)) {
      warning("reducedDims must be a named list. Ignoring.")
      reducedDims <- NULL
    } else {
      reducedDims <- S4Vectors::ReducedDimList(
        lapply(reducedDims, function(x) {
          if (!is.matrix(x) && !is.array(x)) x <- as.matrix(x)
          x
        })
      )
    }
  }

  # Create TreeSummarizedExperiment
  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = assays_list,
    rowData = rowData,
    colData = colData,
    reducedDims = reducedDims,
    rowTree = rowTree,
    rowLinks = rowLinks
  )

  return(tse)
}


#' Generic conversion function for feature tables
#'
#' Automatically detects the input class and converts to the appropriate format.
#' Supports phyloseq, TreeSummarizedExperiment, SingleCellExperiment, and data.frame inputs.
#'
#' @param x Object to convert (phyloseq, TreeSummarizedExperiment, SingleCellExperiment, or data.frame)
#' @param to Target class: "data.frame", "phyloseq", or "TSE" (TreeSummarizedExperiment)
#' @param ... Additional arguments passed to the specific conversion function
#'
#' @return Converted object in the target format
#'
#' @export
#'
#' @examples
#' # Convert from phyloseq to data.frame
#' # df <- convert_feature_table(phylo_obj, to = "data.frame")
#'
#' # Convert from data.frame to phyloseq
#' # ps <- convert_feature_table(df, to = "phyloseq")
#'
#' # Convert from data.frame to TreeSummarizedExperiment
#' # tse <- convert_feature_table(df, to = "TSE")
convert_feature_table <- function(x, to = c("data.frame", "phyloseq", "TSE"), ...) {
  to <- match.arg(to)

  # Detect input class
  x_class <- class(x)[1]

  # Convert to data.frame
  if (to == "data.frame") {
    if (x_class == "phyloseq") {
      return(from_phyloseq(x, ...))
    } else if (x_class %in% c("TreeSummarizedExperiment", "SingleCellExperiment", "SummarizedExperiment")) {
      return(from_TSE(x, ...))
    } else if (is.data.frame(x) || is.matrix(x)) {
      # Already in correct format
      if (is.matrix(x)) {
        return(as.data.frame(x))
      }
      return(x)
    } else {
      stop("Cannot convert ", x_class, " to data.frame. Supported classes: phyloseq, TreeSummarizedExperiment, SingleCellExperiment, data.frame, matrix")
    }

  # Convert to phyloseq
  } else if (to == "phyloseq") {
    if (!requireNamespace("phyloseq", quietly = TRUE)) {
      stop("phyloseq package is required for conversion to phyloseq format")
    }

    if (x_class == "phyloseq") {
      return(x)  # Already phyloseq
    } else if (is.data.frame(x) || is.matrix(x)) {
      return(to_phyloseq(x, ...))
    } else if (x_class %in% c("TreeSummarizedExperiment", "SingleCellExperiment")) {
      # First convert to data.frame, then to phyloseq
      df <- from_TSE(x, ...)
      return(to_phyloseq(df, ...))
    } else {
      stop("Cannot convert ", x_class, " to phyloseq. Supported classes: phyloseq, TreeSummarizedExperiment, SingleCellExperiment, data.frame, matrix")
    }

  # Convert to TreeSummarizedExperiment
  } else if (to == "TSE") {
    if (!requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
      stop("TreeSummarizedExperiment package is required for conversion to TSE format")
    }

    if (x_class %in% c("TreeSummarizedExperiment", "SingleCellExperiment")) {
      return(x)  # Already TSE-like
    } else if (x_class == "phyloseq") {
      # First convert to data.frame, then to TSE
      df <- from_phyloseq(x, ...)
      return(to_TSE(df, ...))
    } else if (is.data.frame(x) || is.matrix(x)) {
      return(to_TSE(x, ...))
    } else {
      stop("Cannot convert ", x_class, " to TreeSummarizedExperiment. Supported classes: TreeSummarizedExperiment, SingleCellExperiment, phyloseq, data.frame, matrix")
    }
  }
}
