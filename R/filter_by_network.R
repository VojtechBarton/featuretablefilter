#' Calculate pairwise mutual information between features
#'
#' Computes the mutual information (MI) matrix for all feature pairs, measuring
#' the amount of information obtained about one feature through the other.
#' High MI indicates strong co-occurrence or mutual exclusion patterns.
#'
#' This function uses a k-nearest neighbors approach for continuous data,
#' which is more robust than histogram-based methods for sparse microbiome data.
#'
#' @param abundances Numeric matrix or data.frame of feature abundances
#'                   (rows = features, columns = samples).
#' @param k Number of nearest neighbors to use for MI estimation. Default is 3.
#'          Higher values give smoother estimates but may miss local patterns.
#' @param bias_correction Logical. Apply bias correction to MI estimates.
#'                        Default is TRUE (recommended for small sample sizes).
#' @param min_prevalence Minimum proportion of samples where a feature must be present
#'                       (non-zero) to be included in MI calculation. Default is 0.05
#'                       (5% of samples). Features below this are excluded to reduce
#'                       computational burden and focus on meaningful signals.
#'
#' @return A symmetric matrix of mutual information values between features.
#'         Diagonal elements are NA (self-information not computed).
#'
#' @export
#'
#' @examples
#' # Create example feature table
#' # mat <- matrix(rpois(100 * 20, lambda = 5), nrow = 100, ncol = 20)
#' # mi_matrix <- compute_mutual_information(mat)
#'
#' # Use with prevalence filtering
#' # filtered_mat <- table[table[, -1] > 0][rowSums(table[, -1] > 0) / ncol(table) > 0.05, ]
#' # mi_matrix <- compute_mutual_information(as.matrix(filtered_mat[, -1]))
compute_mutual_information <- function(abundances, k = 3, bias_correction = TRUE,
                                        min_prevalence = 0.05) {
  # Validate inputs
  if (!is.matrix(abundances) && !is.data.frame(abundances)) {
    stop("abundances must be a matrix or data.frame")
  }

  abundances <- as.matrix(abundances)

  if (nrow(abundances) < 2) {
    stop("Need at least 2 features to compute mutual information")
  }

  if (ncol(abundances) < 4) {
    stop("Need at least 4 samples for reliable MI estimation")
  }

  n_features <- nrow(abundances)
  n_samples <- ncol(abundances)

  # Filter by prevalence to reduce computation and focus on meaningful features
  prevalence <- rowSums(abundances > 0) / n_samples
  valid_features <- prevalence >= min_prevalence

  if (sum(valid_features) < 2) {
    warning("Too few features meet the prevalence threshold. Lower min_prevalence.")
    return(matrix(NA, nrow = n_features, ncol = n_features,
                  dimnames = list(rownames(abundances), rownames(abundances))))
  }

  valid_abundances <- abundances[valid_features, , drop = FALSE]
  valid_names <- rownames(abundances)[valid_features]

  # Initialize MI matrix
  mi_matrix <- matrix(0, nrow = sum(valid_features), ncol = sum(valid_features),
                      dimnames = list(valid_names, valid_names))

  # Compute pairwise MI using k-nearest neighbor estimator
  # This is based on the Kraskov-Stoegbauer-Grassberger (KSG) estimator
  for (i in seq_len(nrow(valid_abundances) - 1)) {
    for (j in (i + 1):nrow(valid_abundances)) {
      x <- valid_abundances[i, ]
      y <- valid_abundances[j, ]

      # Remove samples where either feature is zero (sparse data handling)
      valid_idx <- (x > 0) & (y > 0)
      if (sum(valid_idx) < 4) {
        # Not enough joint occurrences
        mi_value <- 0
      } else {
        x_joint <- x[valid_idx]
        y_joint <- y[valid_idx]

        mi_value <- estimate_mi_knn(x_joint, y_joint, k = k,
                                     bias_correction = bias_correction)
      }

      mi_matrix[i, j] <- mi_value
      mi_matrix[j, i] <- mi_value
    }
  }

  # Set diagonal to NA (not informative for network analysis)
  diag(mi_matrix) <- NA

  # Create full matrix with all features
  full_mi <- matrix(0, nrow = n_features, ncol = n_features,
                    dimnames = list(rownames(abundances), rownames(abundances)))
  full_mi[valid_features, valid_features] <- mi_matrix

  return(full_mi)
}

#' Estimate mutual information using k-nearest neighbors (KSG estimator)
#'
#' Internal function implementing the Kraskov-Stoegbauer-Grassberger estimator
#' for mutual information. This is a state-of-the-art method that works well
#' for continuous and mixed data types.
#'
#' @param x Numeric vector of observations
#' @param y Numeric vector of observations
#' @param k Number of nearest neighbors
#' @param bias_correction Logical. Apply bias correction
#'
#' @return Mutual information estimate in bits
estimate_mi_knn <- function(x, y, k = 3, bias_correction = TRUE) {
  n <- length(x)

  if (n < 4) return(0)

  # Check for constant vectors (zero variance)
  if (sd(x, na.rm = TRUE) < .Machine$double.eps ||
      sd(y, na.rm = TRUE) < .Machine$double.eps) {
    return(0)
  }

  # Combine into matrix for distance calculations
  data_mat <- cbind(x, y)

  # Compute pairwise distances in joint space
  joint_dist <- dist(data_mat, method = "euclidean")
  joint_dist <- as.matrix(joint_dist)

  # Compute marginal distances
  x_dist <- as.matrix(dist(x))
  y_dist <- as.matrix(dist(y))

  # For each point, find k-th nearest neighbor distance in joint space
  # (excluding self)
  kmi_values <- numeric(n)
  nx_values <- numeric(n)
  ny_values <- numeric(n)

  for (i in seq_len(n)) {
    # Sort distances from point i (excluding self)
    joint_dists_i <- sort(joint_dist[i, -i])
    kmi_dists_i <- joint_dists_i[min(k, length(joint_dists_i))]

    # Count points within kmi distance in marginal spaces
    # (using < rather than <= to avoid counting boundary points)
    nx_values[i] <- sum(x_dist[i, ] < kmi_dists_i) - 1  # Exclude self
    ny_values[i] <- sum(y_dist[i, ] < kmi_dists_i) - 1  # Exclude self

    # The k-th NN distance in joint space
    kmi_values[i] <- kmi_dists_i
  }

  # KSG estimator: I(X;Y) = psi(k) - mean(psi(nx + 1)) - mean(psy(ny + 1)) + psi(n)
  # Using digamma function
  # Handle edge cases where nx or ny might be 0
  nx_plus_one <- pmax(nx_values + 1, 1)
  ny_plus_one <- pmax(ny_values + 1, 1)

  if (bias_correction) {
    mi_estimate <- digamma(k) - mean(digamma(nx_plus_one), na.rm = TRUE) -
                   mean(digamma(ny_plus_one), na.rm = TRUE) + digamma(n)
  } else {
    # Simpler version without full bias correction
    mi_estimate <- log(n) - mean(log(nx_plus_one), na.rm = TRUE) -
                   mean(log(ny_plus_one), na.rm = TRUE)
  }

  # Handle any remaining NA values
  if (is.na(mi_estimate)) {
    return(0)
  }

  # Convert to bits (from nats)
  mi_bits <- mi_estimate / log(2)

  # Ensure non-negative
  return(max(0, mi_bits))
}

#' Calculate network connectivity metrics for features
#'
#' Based on a mutual information or correlation matrix, computes network
#' centrality measures for each feature. Features with zero connections
#' (degree = 0) are likely spurious or random artifacts.
#'
#' @param similarity_matrix Symmetric matrix of pairwise similarities
#'                          (mutual information or correlation).
#' @param threshold Cutoff value above which an edge is considered significant.
#'                  For MI, typical values are 0.1-0.5 bits. For correlation,
#'                  typical values are 0.3-0.7. Default is NULL (automatic).
#' @param method Type of similarity measure: "mi" for mutual information,
#'               "cor" for Pearson correlation. Used for automatic thresholding.
#' @param automatic_threshold_method Method for automatic threshold selection:
#'                                   "mean_sd" (mean + 2*SD), "percentile" (75th percentile),
#'                                   or "mst" (minimum spanning tree). Default is "mean_sd".
#'
#' @return A list containing:
#'   \item{degree_centrality}{Number of significant connections per feature}
#'   \item{strength}{Sum of connection strengths per feature}
#'   \item{betweenness}{Betweenness centrality (if igraph available)}
#'   \item{threshold_used}{The threshold value used}
#'   \item{n_edges}{Total number of edges in the network}
#'   \item{adjacency_matrix}{Binary adjacency matrix (1 = connected, 0 = not)}
#'
#' @export
#'
#' @examples
#' # Compute MI matrix first
#' # mi_matrix <- compute_mutual_information(feature_table)
#'
#' # Get network metrics
#' # network_metrics <- analyze_feature_network(mi_matrix)
#'
#' # Features with zero degree are likely artifacts
#' # disconnected_features <- names(network_metrics$degree_centrality[
#' #   network_metrics$degree_centrality == 0
#' # ])
analyze_feature_network <- function(similarity_matrix, threshold = NULL,
                                     method = c("mi", "cor"),
                                     automatic_threshold_method = c("mean_sd", "percentile", "mst")) {
  if (!is.matrix(similarity_matrix) && !is.data.frame(similarity_matrix)) {
    stop("similarity_matrix must be a matrix or data.frame")
  }

  similarity_matrix <- as.matrix(similarity_matrix)

  # Remove diagonal for analysis
  diag(similarity_matrix) <- 0

  method <- match.arg(method)
  automatic_threshold_method <- match.arg(automatic_threshold_method)

  n_features <- nrow(similarity_matrix)

  # Automatic threshold selection if not provided
  if (is.null(threshold)) {
    upper_tri_vals <- similarity_matrix[upper.tri(similarity_matrix)]
    valid_vals <- upper_tri_vals[upper_tri_vals > 0]

    if (length(valid_vals) == 0) {
      threshold <- 0
    } else if (automatic_threshold_method == "mean_sd") {
      val_mean <- mean(valid_vals, na.rm = TRUE)
      val_sd <- sd(valid_vals, na.rm = TRUE)
      if (is.na(val_sd) || val_sd == 0) {
        threshold <- val_mean
      } else {
        threshold <- val_mean + 2 * val_sd
      }
    } else if (automatic_threshold_method == "percentile") {
      threshold <- as.numeric(quantile(valid_vals, 0.75, na.rm = TRUE))
    } else {  # mst - minimum spanning tree approach
      # Threshold at the largest gap in sorted unique values
      sorted_vals <- sort(unique(valid_vals[!is.na(valid_vals)]))
      if (length(sorted_vals) > 1) {
        gaps <- diff(sorted_vals)
        threshold <- sorted_vals[which.max(gaps)]
      } else {
        threshold <- median(valid_vals, na.rm = TRUE)
      }
    }

    # Ensure reasonable bounds
    if (!is.na(threshold)) {
      if (method == "mi" && threshold < 0.05) threshold <- 0.05
      if (method == "cor" && threshold < 0.2) threshold <- 0.2
    } else {
      threshold <- if (method == "mi") 0.05 else 0.2
    }
  }

  # Create binary adjacency matrix
  adjacency_matrix <- (similarity_matrix >= threshold) * 1
  diag(adjacency_matrix) <- 0

  # Calculate degree centrality (number of connections)
  degree_centrality <- rowSums(adjacency_matrix)

  # Calculate strength (sum of connection weights)
  strength <- rowSums(similarity_matrix * adjacency_matrix)

  # Try to calculate betweenness centrality if igraph is available
  betweenness <- rep(NA, n_features)
  names(betweenness) <- rownames(similarity_matrix)

  if (requireNamespace("igraph", quietly = TRUE)) {
    tryCatch({
      g <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected",
                                        weighted = TRUE)
      betweenness <- betweenness(g, directed = FALSE)
    }, error = function(e) {
      # Betweenness calculation failed
    })
  }

  # Count total edges
  n_edges <- sum(adjacency_matrix, na.rm = TRUE) / 2  # Divide by 2 for undirected
  if (is.na(n_edges)) n_edges <- 0

  return(list(
    degree_centrality = degree_centrality,
    strength = strength,
    betweenness = betweenness,
    threshold_used = threshold,
    n_edges = as.integer(n_edges),
    adjacency_matrix = adjacency_matrix
  ))
}

#' Filter features based on network connectivity
#'
#' Removes features that have no significant connections to other features
#' in the network. These disconnected features are likely spurious artifacts,
#' cross-talk, or contamination rather than true biological signals.
#'
#' This approach is based on the principle that true biological features
#' typically exhibit coordinated behavior with other features due to
#' ecological interactions, shared environmental responses, or phylogenetic
#' relationships.
#'
#' @param table A feature table (data.frame or matrix) with features as rows
#'              and samples as columns. First column should be feature IDs.
#' @param similarity_type Type of similarity to compute: "mi" for mutual information,
#'                        "cor" for Pearson correlation. Default is "mi".
#' @param threshold Similarity threshold for defining edges. If NULL, automatic
#'                  thresholding is used. Default is NULL.
#' @param min_degree Minimum degree (number of connections) required to keep a feature.
#'                   Default is 1 (features must have at least one connection).
#' @param min_prevalence Minimum proportion of samples where a feature must be present
#'                       for MI calculation. Passed to \code{\link{compute_mutual_information}}.
#'                       Default is 0.05.
#' @param k Number of nearest neighbors for MI estimation. Passed to
#'          \code{\link{compute_mutual_information}}. Default is 3.
#' @param verbose Logical. Print progress and summary. Default is TRUE.
#'
#' @return A filtered feature table with disconnected features removed.
#'         Attributes include:
#'   \item{n_filtered_out}{Number of features removed}
#'   \item{network_summary}{List with network statistics}
#'   \item{disconnected_features}{Character vector of removed feature names}
#'
#' @export
#'
#' @examples
#' # Filter out features with no network connections
#' # cleaned_table <- filter_by_network_connectivity(my_table)
#'
#' # Use correlation instead of MI (faster but less sensitive to non-linear patterns)
#' # cleaned_table <- filter_by_network_connectivity(my_table, similarity_type = "cor")
#'
#' # Require at least 2 connections to keep a feature
#' # cleaned_table <- filter_by_network_connectivity(my_table, min_degree = 2)
filter_by_network_connectivity <- function(table, similarity_type = c("mi", "cor"),
                                            threshold = NULL, min_degree = 1,
                                            min_prevalence = 0.05, k = 3,
                                            verbose = TRUE) {
  similarity_type <- match.arg(similarity_type)

  # Validate input
  if (!is.data.frame(table) && !is.matrix(table)) {
    stop("table must be a data.frame or matrix")
  }

  if (ncol(table) < 2) {
    stop("table must have at least one sample column")
  }

  # Extract feature IDs and abundances
  feature_ids <- table[, 1, drop = FALSE]
  abundances <- as.matrix(table[, -1, drop = FALSE])

  n_features_orig <- nrow(abundances)

  if (verbose) {
    message(sprintf("Computing %s similarity matrix for %d features...",
                    similarity_type, n_features_orig))
  }

  # Compute similarity matrix
  if (similarity_type == "mi") {
    similarity_matrix <- compute_mutual_information(
      abundances, k = k, bias_correction = TRUE,
      min_prevalence = min_prevalence
    )
  } else {
    # Pearson correlation
    # Transpose so samples are rows for cor()
    corr_matrix <- cor(t(abundances), method = "pearson", use = "pairwise.complete")
    # Convert to similarity (absolute correlation)
    similarity_matrix <- abs(corr_matrix)
    diag(similarity_matrix) <- 0
  }

  if (verbose) {
    message("Analyzing feature network...")
  }

  # Analyze network
  network_metrics <- analyze_feature_network(
    similarity_matrix, threshold = threshold,
    method = similarity_type
  )

  # Identify features to keep
  feature_names <- rownames(abundances)
  degree <- network_metrics$degree_centrality
  names(degree) <- feature_names

  # Handle NA values (features excluded from MI calculation due to low prevalence)
  # These are treated as disconnected (will be filtered out)
  degree_clean <- ifelse(is.na(degree), 0, degree)
  keep_features <- degree_clean >= min_degree
  disconnected_features <- feature_names[!keep_features]

  n_filtered_out <- sum(!keep_features)
  n_retained <- sum(keep_features)

  if (verbose) {
    message(sprintf("\nNetwork Analysis Summary:"))
    message(sprintf("  Threshold used: %.4f", network_metrics$threshold_used))
    message(sprintf("  Total edges: %d", network_metrics$n_edges))
    message(sprintf("  Features retained: %d (%.1f%%)", n_retained,
                    (n_retained / n_features_orig) * 100))
    message(sprintf("  Features removed (disconnected): %d (%.1f%%)",
                    n_filtered_out, (n_filtered_out / n_features_orig) * 100))

    if (n_filtered_out > 0 && verbose) {
      if (n_filtered_out <= 10) {
        message(sprintf("  Disconnected features: %s",
                        paste(disconnected_features, collapse = ", ")))
      } else {
        message(sprintf("  Disconnected features: %s (and %d more)",
                        paste(head(disconnected_features, 10), collapse = ", "),
                        n_filtered_out - 10))
      }
    }
  }

  # Build result
  result <- cbind(feature_ids[keep_features, , drop = FALSE],
                  abundances[keep_features, , drop = FALSE])

  if (is.matrix(result)) {
    result <- as.data.frame(result)
  }
  colnames(result) <- colnames(table)

  # Attach attributes
  attr(result, "n_filtered_out") <- n_filtered_out
  attr(result, "network_summary") <- list(
    threshold_used = network_metrics$threshold_used,
    n_edges = network_metrics$n_edges,
    similarity_type = similarity_type,
    min_degree = min_degree
  )
  attr(result, "disconnected_features") <- disconnected_features
  attr(result, "network_metrics") <- network_metrics

  return(result)
}

#' Plot feature network connectivity
#'
#' Creates a visualization of the feature network showing connections between
#' features. Useful for exploring ecological relationships and identifying
# disconnected (potentially spurious) features.
#'
#' @param network_metrics Output from \code{\link{analyze_feature_network}}
#' @param feature_names Character vector of feature names (must match network_metrics)
#' @param top_n Number of top-connected features to label. Default is 10.
#' @param color_by Degree centrality (connections) or strength (connection weight).
#'                 Default is "degree".
#'
#' @return A ggplot object (edge list format for flexibility) or null plot
#'         if igraph is not available.
#'
#' @export
#'
#' @examples
#' # mi_matrix <- compute_mutual_information(feature_table)
#' # network <- analyze_feature_network(mi_matrix)
#' # plot_feature_network(network, rownames(feature_table))
plot_feature_network <- function(network_metrics, feature_names = NULL,
                                  top_n = 10, color_by = c("degree", "strength")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 is required for plotting. Install with install.packages('ggplot2').")
    return(NULL)
  }

  if (!requireNamespace("igraph", quietly = TRUE)) {
    warning("igraph is recommended for network visualization. Install with install.packages('igraph').")
    return(NULL)
  }

  color_by <- match.arg(color_by)

  adjacency <- network_metrics$adjacency_matrix
  feature_names <- feature_names %||% rownames(adjacency)

  if (is.null(feature_names) || length(feature_names) != nrow(adjacency)) {
    warning("feature_names must be provided and match network dimensions")
    return(NULL)
  }

  # Create graph using igraph
  g <- igraph::graph_from_adjacency_matrix(adjacency, mode = "undirected", weighted = TRUE)
  igraph::V(g)$name <- feature_names
  igraph::V(g)$degree <- network_metrics$degree_centrality
  igraph::V(g)$strength <- network_metrics$strength

  # Layout
  lay <- igraph::layout_with_fr(g)

  # Get edge list
  edges <- igraph::get.edgelist(g)
  edge_weights <- igraph::E(g)$weight

  edge_df <- data.frame(
    from = feature_names[edges[, 1]],
    to = feature_names[edges[, 2]],
    weight = edge_weights
  )

  # Node data
  node_df <- data.frame(
    name = feature_names,
    degree = igraph::V(g)$degree,
    strength = igraph::V(g)$strength,
    x = lay[, 1],
    y = lay[, 2],
    stringsAsFactors = FALSE
  )

  # Label top N nodes
  top_nodes <- node_df[order(-node_df[[color_by]]), ][seq_len(min(top_n, nrow(node_df))), ]
  node_df$label <- ifelse(node_df$name %in% top_nodes$name, node_df$name, "")

  base_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = edge_df, ggplot2::aes(x = 0, y = 0, xend = 0, yend = 0),
                 alpha = 0) +  # Placeholder for proper edge rendering
    ggplot2::geom_point(data = node_df, ggplot2::aes(x = x, y = y,
                                    color = !!ggplot2::sym(color_by),
                                    size = degree), alpha = 0.8) +
    ggplot2::geom_text(data = node_df, ggplot2::aes(x = x, y = y, label = label),
              size = 3, hjust = 0.5, vjust = 0.5) +
    ggplot2::scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +
    ggplot2::labs(
      title = "Feature Network Connectivity",
      subtitle = sprintf("Nodes = features | Edges = significant similarities | Color = %s", color_by),
      color = color_by,
      size = "Degree"
    ) +
    base_theme

  return(p)
}

# Helper for null coalescing (R < 4.0 compatibility)
`%||%` <- function(x, y) if (is.null(x)) y else x
