#' Compute quality control metrics comparing original and filtered tables
#'
#' Calculates multiple QC metrics to assess the impact of filtering on a feature table.
#' Includes sparsity change, retention rates, rank-abundance stability, and compositional similarity.
#'
#' For Rank-Abundance Stability, top N features are independently identified in each table by:
#' (1) summing feature abundances across all samples, (2) converting to relative abundances,
#' (3) selecting top N, then (4) comparing overlap and rank correlation of common features.
#'
#' @param original_table Original feature table (data.frame) before filtering.
#' @param filtered_table Filtered feature table (data.frame) after filtering.
#' @param top_n Number of taxa to consider for Rank-Abundance Stability. Default is 10.
#'
#' @return A list containing:
#'   \item{sparsity_original}{Sparsity of original table}
#'   \item{sparsity_filtered}{Sparsity of filtered table}
#'   \item{sparsity_drop_percent}{Percentage point drop in sparsity}
#'   \item{read_retention_percent}{Percentage of reads retained}
#'   \item{feature_retention_percent}{Percentage of features retained}
#'   \item{sample_retention_percent}{Percentage of samples retained}
#'   \item{rank_abundance_correlation}{Spearman correlation of ranks for features in both top N lists}
#'   \item{rank_abundance_pvalue}{P-value for rank-abundance correlation}
#'   \item{top_n_overlap_count}{Number of features appearing in both top N lists}
#'   \item{top_n_overlap_percent}{Percentage of top N features that overlap}
#'   \item{top_n_jaccard_similarity}{Jaccard similarity index of top N feature sets}
#'   \item{orig_top_features}{Character vector of top N features in original table}
#'   \item{orig_top_rels}{Relative abundances of top N features in original table}
#'   \item{filt_top_features}{Character vector of top N features in filtered table}
#'   \item{filt_top_rels}{Relative abundances of top N features in filtered table}
#'   \item{procrustes_m2}{Procrustes M^2 statistic (Bray-Curtis based)}
#'   \item{procrustes_correlation}{Procrustes correlation coefficient}
#'
#' @export
#'
#' @examples
#' # Compare original and filtered tables
#' # qc <- compute_filtering_qc(original_table, filtered_table)
#' # print(qc)
compute_filtering_qc <- function(original_table, filtered_table, top_n = 10) {
  # Extract abundance matrices
  orig_abund <- as.matrix(original_table[, -1, drop = FALSE])
  filt_abund <- as.matrix(filtered_table[, -1, drop = FALSE])
  orig_features <- original_table[, 1]
  filt_features <- filtered_table[, 1]
  orig_samples <- colnames(orig_abund)
  filt_samples <- colnames(filt_abund)

  # Find common features and samples for comparison
  common_features <- intersect(orig_features, filt_features)
  common_samples <- intersect(orig_samples, filt_samples)

  # === Sparsity ===
  sparsity_orig <- sum(orig_abund == 0) / length(orig_abund)
  sparsity_filt <- sum(filt_abund == 0) / length(filt_abund)
  sparsity_drop <- (sparsity_orig - sparsity_filt) * 100

  # === Retention rates ===
  total_reads_orig <- sum(orig_abund)
  total_reads_filt <- sum(filt_abund)
  read_retention <- (total_reads_filt / total_reads_orig) * 100

  feature_retention <- (length(common_features) / length(orig_features)) * 100
  sample_retention <- (length(common_samples) / length(orig_samples)) * 100

  # === Rank-Abundance Stability (Top N Taxa Invariance) ===
  # Correct approach: independently find top N in each table, then compare

  # Calculate total abundance per feature (sum across all samples)
  orig_total <- rowSums(orig_abund)
  filt_total <- rowSums(filt_abund)

  # Convert to relative abundances (within each table)
  orig_rel_total <- orig_total / sum(orig_total)
  filt_rel_total <- filt_total / sum(filt_total)

  # Get top N features independently from each table
  orig_top_idx <- order(orig_rel_total, decreasing = TRUE)[seq_len(min(top_n, length(orig_rel_total)))]
  filt_top_idx <- order(filt_rel_total, decreasing = TRUE)[seq_len(min(top_n, length(filt_rel_total)))]

  orig_top_features <- orig_features[orig_top_idx]
  filt_top_features <- filt_features[filt_top_idx]

  # Store relative abundances of top features for reporting
  orig_top_rels <- orig_rel_total[orig_top_idx]
  filt_top_rels <- filt_rel_total[filt_top_idx]

  # Calculate overlap metrics
  overlap_count <- length(intersect(orig_top_features, filt_top_features))
  jaccard_similarity <- overlap_count / length(union(orig_top_features, filt_top_features))
  overlap_percent <- (overlap_count / top_n) * 100

  # Spearman correlation: compare ranks of common features in top N lists
  common_in_top <- intersect(orig_top_features, filt_top_features)
  if (length(common_in_top) >= 2) {
    # Get ranks of common features in original top N
    orig_ranks <- sapply(common_in_top, function(f) {
      which(orig_top_features == f)
    })
    # Get ranks of common features in filtered top N
    filt_ranks <- sapply(common_in_top, function(f) {
      which(filt_top_features == f)
    })

    rank_corr_test <- cor.test(as.numeric(orig_ranks), as.numeric(filt_ranks), method = "spearman")
    rank_abundance_correlation <- rank_corr_test$estimate
    rank_abundance_pvalue <- rank_corr_test$p.value
  } else {
    rank_abundance_correlation <- NA
    rank_abundance_pvalue <- NA
  }

  # === Procrustes Analysis (Bray-Curtis based) ===
  procrustes_m2 <- NA
  procrustes_correlation <- NA

  if (requireNamespace("vegan", quietly = TRUE)) {
    tryCatch({
      # Create Bray-Curtis distance matrices
      # Subset original and filtered tables to common features and samples
      orig_subset <- orig_abund[match(common_features, orig_features),
                                 match(common_samples, orig_samples), drop = FALSE]
      filt_subset <- filt_abund[match(common_features, filt_features),
                                 match(common_samples, filt_samples), drop = FALSE]

      veg_dist_orig <- vegan::vegdist(t(orig_subset), method = "bray")
      veg_dist_filt <- vegan::vegdist(t(filt_subset), method = "bray")

      # Perform Procrustes analysis on distance matrices
      proc <- vegan::procrustes(veg_dist_orig, veg_dist_filt, scale = TRUE)

      # Calculate M² and r² from ss (sum of squared residuals)
      # M² = ss / sum(Y²), r² = 1 - M²
      if (!is.null(proc$ss)) {
        y_sum_sq <- sum(proc$Yrot^2)
        if (y_sum_sq > 0) {
          procrustes_m2 <- proc$ss / y_sum_sq
          procrustes_correlation <- 1 - procrustes_m2
        }
      }

    }, error = function(e) {
      warning("Procrustes analysis failed: ", e$message)
    })
  } else {
    warning("Package 'vegan' required for Procrustes analysis. Install with: install.packages('vegan')")
  }

  return(list(
    sparsity_original = sparsity_orig,
    sparsity_filtered = sparsity_filt,
    sparsity_drop_percent = sparsity_drop,
    read_retention_percent = read_retention,
    feature_retention_percent = feature_retention,
    sample_retention_percent = sample_retention,
    rank_abundance_correlation = rank_abundance_correlation,
    rank_abundance_pvalue = rank_abundance_pvalue,
    top_n_overlap_count = overlap_count,
    top_n_overlap_percent = overlap_percent,
    top_n_jaccard_similarity = jaccard_similarity,
    orig_top_features = orig_top_features,
    orig_top_rels = orig_top_rels,
    filt_top_features = filt_top_features,
    filt_top_rels = filt_top_rels,
    procrustes_m2 = procrustes_m2,
    procrustes_correlation = procrustes_correlation
  ))
}
