#' Compute quality control metrics comparing original and filtered tables
#'
#' Calculates multiple QC metrics to assess the impact of filtering on a feature table.
#' Includes sparsity change, retention rates, rank-abundance stability, compositional similarity,
#' and diversity retention (Effective Number of Species / Hill numbers).
#'
#' For Rank-Abundance Stability, top N features are independently identified in each table by:
#' (1) summing feature abundances across all samples, (2) converting to relative abundances,
#' (3) selecting top N, then (4) comparing overlap and rank correlation of common features.
#'
#' For Effective Number of Species (ENS), computes Hill numbers (Shannon and Simpson diversity
#' converted to effective species count) per sample and reports how much diversity is retained
#' after filtering. This reveals whether filtering preserves biological diversity or flattens it.
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
#'   \item{shannon_ens_original}{Mean Shannon effective number of species (Hill q=1) across samples}
#'   \item{shannon_ens_filtered}{Mean Shannon ENS after filtering}
#'   \item{shannon_ens_retention_percent}{Percentage of Shannon ENS retained}
#'   \item{simpson_ens_original}{Mean Simpson effective number of species (Hill q=2) across samples}
#'   \item{simpson_ens_filtered}{Mean Simpson ENS after filtering}
#'   \item{simpson_ens_retention_percent}{Percentage of Simpson ENS retained}
#'
#' @export
#'
#' @examples
#' # Compare original and filtered tables
#' # qc <- compute_filtering_qc(original_table, filtered_table)
#' # print(qc)
#'
#' # Access diversity retention metrics
#' # qc$shannon_ens_retention_percent  # How much Shannon diversity (q=1) was retained
#' # qc$simpson_ens_retention_percent  # How much Simpson diversity (q=2) was retained
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

  # === Effective Number of Species (Hill Numbers) ===
  # Shannon diversity (q=1): exp(H') where H' is Shannon entropy
  # Simpson diversity (q=2): 1 / (1 - D) where D is Simpson's index
  # These convert diversity indices to "effective number of species"

  # Calculate Shannon ENS per sample (exp(Shannon entropy))
  calc_shannon_ens <- function(abund_matrix) {
    apply(abund_matrix, 2, function(sample_vec) {
      counts <- sample_vec[sample_vec > 0]
      if (length(counts) <= 1) return(NA)
      total <- sum(counts)
      probs <- counts / total
      shannon_entropy <- -sum(probs * log(probs))
      exp(shannon_entropy)  # Convert to effective number of species
    })
  }

  # Calculate Simpson ENS per sample (1 / sum(p^2))
  # This is the inverse Simpson index, equivalent to Hill number q=2
  calc_simpson_ens <- function(abund_matrix) {
    apply(abund_matrix, 2, function(sample_vec) {
      counts <- sample_vec[sample_vec > 0]
      if (length(counts) <= 1) return(NA)
      total <- sum(counts)
      probs <- counts / total
      simpson_d <- sum(probs^2)
      if (simpson_d >= 1) return(1)  # Only one species present
      1 / simpson_d  # Inverse Simpson = effective number of species
    })
  }

  # Compute ENS metrics for common samples only
  orig_common_samples <- orig_abund[, match(common_samples, orig_samples), drop = FALSE]
  filt_common_samples <- filt_abund[, match(common_samples, filt_samples), drop = FALSE]

  # Shannon ENS
  shannon_ens_orig <- calc_shannon_ens(orig_common_samples)
  shannon_ens_filt <- calc_shannon_ens(filt_common_samples)

  # Simpson ENS
  simpson_ens_orig <- calc_simpson_ens(orig_common_samples)
  simpson_ens_filt <- calc_simpson_ens(filt_common_samples)

  # Aggregate: mean ENS across samples (excluding NAs)
  shannon_ens_original <- mean(shannon_ens_orig, na.rm = TRUE)
  shannon_ens_filtered <- mean(shannon_ens_filt, na.rm = TRUE)
  simpson_ens_original <- mean(simpson_ens_orig, na.rm = TRUE)
  simpson_ens_filtered <- mean(simpson_ens_filt, na.rm = TRUE)

  # Calculate retention percentages
  if (!is.na(shannon_ens_original) && shannon_ens_original > 0) {
    shannon_ens_retention_percent <- (shannon_ens_filtered / shannon_ens_original) * 100
  } else {
    shannon_ens_retention_percent <- NA
  }

  if (!is.na(simpson_ens_original) && simpson_ens_original > 0) {
    simpson_ens_retention_percent <- (simpson_ens_filtered / simpson_ens_original) * 100
  } else {
    simpson_ens_retention_percent <- NA
  }

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

  # Spearman and Pearson correlation: compare abundances of common features in top N lists
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

    # Spearman correlation (rank-based)
    rank_corr_test <- cor.test(as.numeric(orig_ranks), as.numeric(filt_ranks), method = "spearman")
    rank_abundance_correlation <- rank_corr_test$estimate
    rank_abundance_pvalue <- rank_corr_test$p.value

    # Pearson correlation (linear relationship of relative abundances)
    orig_abunds <- sapply(common_in_top, function(f) {
      orig_top_rels[which(orig_top_features == f)]
    })
    filt_abunds <- sapply(common_in_top, function(f) {
      filt_top_rels[which(filt_top_features == f)]
    })

    pearson_corr_result <- tryCatch({
      ct <- cor.test(as.numeric(orig_abunds), as.numeric(filt_abunds), method = "pearson")
      list(estimate = ct$estimate, p.value = ct$p.value)
    }, error = function(e) {
      list(estimate = NA_real_, p.value = NA_real_)
    })
    pearson_abundance_correlation <- pearson_corr_result$estimate
    pearson_abundance_pvalue <- pearson_corr_result$p.value
  } else {
    rank_abundance_correlation <- NA
    rank_abundance_pvalue <- NA
    pearson_abundance_correlation <- NA
    pearson_abundance_pvalue <- NA
  }

  # === Procrustes Analysis (Bray-Curtis based) ===
  procrustes_m2 <- NA
  procrustes_correlation <- NA
  procrustes_pvalue <- NA

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

      # Perform permutation test for p-value using protest
      # This tests the significance of the correlation between two configurations
      protest_result <- vegan::protest(veg_dist_orig, veg_dist_filt, permutations = 999)
      if (!is.null(protest_result$signif)) {
        procrustes_pvalue <- protest_result$signif
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
    pearson_abundance_correlation = pearson_abundance_correlation,
    pearson_abundance_pvalue = pearson_abundance_pvalue,
    top_n_overlap_count = overlap_count,
    top_n_overlap_percent = overlap_percent,
    top_n_jaccard_similarity = jaccard_similarity,
    orig_top_features = orig_top_features,
    orig_top_rels = orig_top_rels,
    filt_top_features = filt_top_features,
    filt_top_rels = filt_top_rels,
    procrustes_m2 = procrustes_m2,
    procrustes_correlation = procrustes_correlation,
    procrustes_pvalue = procrustes_pvalue,
    shannon_ens_original = shannon_ens_original,
    shannon_ens_filtered = shannon_ens_filtered,
    shannon_ens_retention_percent = shannon_ens_retention_percent,
    simpson_ens_original = simpson_ens_original,
    simpson_ens_filtered = simpson_ens_filtered,
    simpson_ens_retention_percent = simpson_ens_retention_percent
  ))
}

#' Calculate Shannon Effective Number of Species (Hill number q=1)
#'
#' Converts Shannon entropy to effective number of species (True Diversity).
#' For a community with perfect evenness among S species, ENS = S.
#' For uneven communities, ENS < S, reflecting that rare species contribute less
#' to the "effective" diversity.
#'
#' @param abund_matrix Numeric matrix of abundances (rows = features, columns = samples).
#' @return Numeric vector of Shannon ENS values, one per sample.
#' @export
#'
#' @examples
#' # Example feature table (rows = samples, cols = features)
#' mat <- matrix(c(25, 25, 25, 25, 80, 5, 5, 5), nrow = 2, byrow = TRUE)
#' colnames(mat) <- c("f1", "f2", "f3", "f4")
#' calc_shannon_ens(mat)
calc_shannon_ens <- function(abund_matrix) {
  apply(abund_matrix, 2, function(sample_vec) {
    counts <- sample_vec[sample_vec > 0]
    if (length(counts) <= 1) return(NA)
    total <- sum(counts)
    probs <- counts / total
    shannon_entropy <- -sum(probs * log(probs))
    exp(shannon_entropy)  # Convert to effective number of species
  })
}

#' Calculate Simpson Effective Number of Species (Hill number q=2)
#'
#' Converts Simpson's index to effective number of species using the inverse
#' Simpson formula (1/sum(p^2)). This metric is more sensitive to dominant
#' species than Shannon ENS.
#'
#' @param abund_matrix Numeric matrix of abundances (rows = features, columns = samples).
#' @return Numeric vector of Simpson ENS values, one per sample.
#' @export
#'
#' @examples
#' # Example feature table
#' mat <- matrix(c(20, 20, 20, 20, 20, 50, 10, 10, 15, 15), nrow = 2, byrow = TRUE)
#' colnames(mat) <- c("f1", "f2", "f3", "f4", "f5")
#' calc_simpson_ens(mat)
calc_simpson_ens <- function(abund_matrix) {
  apply(abund_matrix, 2, function(sample_vec) {
    counts <- sample_vec[sample_vec > 0]
    if (length(counts) == 0) return(NA)  # No data
    if (length(counts) == 1) return(1)   # Single species = ENS of 1
    total <- sum(counts)
    probs <- counts / total
    simpson_d <- sum(probs^2)
    1 / simpson_d  # Inverse Simpson = effective number of species
  })
}

#' Calculate Hill Numbers (Effective Number of Species) for a range of q values
#'
#' Computes diversity profiles across different sensitivity parameters (q).
#' Lower q values are more sensitive to rare species; higher q values focus
#' on dominant species.
#'
#' @param abund_matrix Numeric matrix of abundances (rows = features, columns = samples).
#' @param q Vector of Hill number order parameters. Default includes q=0 (richness),
#'   q=1 (Shannon ENS), and q=2 (Simpson ENS).
#' @return Matrix with rows = samples and columns = q values, containing ENS estimates.
#' @export
#'
#' @examples
#' mat <- matrix(c(25, 25, 25, 25, 80, 5, 5, 5), nrow = 2, byrow = TRUE)
#' colnames(mat) <- c("f1", "f2", "f3", "f4")
#' calc_hill_numbers(mat, q = c(0, 0.5, 1, 2, 3))
calc_hill_numbers <- function(abund_matrix, q = c(0, 1, 2)) {
  result <- sapply(q, function(q_val) {
    apply(abund_matrix, 2, function(sample_vec) {
      counts <- sample_vec[sample_vec > 0]
      if (length(counts) == 0) return(NA)  # No data
      if (length(counts) == 1) {
        # Single species: ENS = 1 for all q > 0, richness = 1 for q = 0
        if (q_val == 0) return(1)
        return(1)
      }
      total <- sum(counts)
      probs <- counts / total

      if (q_val == 0) {
        # q=0: Species richness (count of non-zero features)
        length(counts)
      } else if (q_val == 1) {
        # q=1: Shannon ENS (exp of Shannon entropy)
        shannon_entropy <- -sum(probs * log(probs))
        exp(shannon_entropy)
      } else {
        # q>0, q!=1: Generalized Hill number formula
        # ^qD = (sum(p_i^q))^(1/(1-q))
        sum_p_q <- sum(probs^q_val)
        if (sum_p_q == 0) return(NA)
        sum_p_q^(1 / (1 - q_val))
      }
    })
  }, simplify = "matrix")
  # Ensure result is always a matrix even with single column
  if (!is.matrix(result)) {
    result <- matrix(result, ncol = length(q))
    colnames(result) <- q
  }
  result
}
