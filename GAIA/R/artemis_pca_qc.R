#' PC-Metadata Association Test
#'
#' @description Runs PCA on a log2-scale matrix, then tests the association of
#' each principal component with each metadata variable. Returns a p-value
#' matrix that summarises which variables explain sample-level variance.
#'
#' Intended use: run after loading and filtering, before differential analysis,
#' to distinguish expected biological drivers (group, timepoint) from
#' unexpected technical ones (plate, batch). If a technical variable is
#' significant, consider correcting with \code{POSEIDON} before proceeding.
#'
#' @param matrix Numeric matrix, features x samples (log2-scale). Must have no
#'   NA values — filter with \code{HADES_filter_olink_proteins()} first if needed.
#'   Rownames = feature IDs, colnames = sample IDs.
#' @param sample_meta Data.frame with one row per sample. Matched to
#'   \code{colnames(matrix)} by \code{sample_col} (or rownames if NULL).
#' @param n_pcs Integer. Number of top PCs to test. Default: 10.
#' @param sample_col Character or NULL. Column in \code{sample_meta} holding
#'   sample IDs matching \code{colnames(matrix)}. NULL = use rownames.
#'   Default: NULL.
#' @param metadata_cols Character vector or NULL. Columns of \code{sample_meta}
#'   to test. NULL = all columns except \code{sample_col}. Default: NULL.
#' @param ntop Integer or NULL. Number of most variable features (by row
#'   variance) to use for PCA. NULL = use all features. Default: NULL.
#' @param scale Logical. Scale features to unit variance before PCA. Recommended
#'   for NPX data where proteins have different dynamic ranges. Default: TRUE.
#' @param categorical_threshold Integer. Variables with at most this many unique
#'   non-NA values are treated as categorical (Kruskal-Wallis). Variables with
#'   more unique values are treated as continuous (Spearman). Default: 10.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return An S3 object of class \code{"artemis_pc_assoc"} containing:
#' \describe{
#'   \item{pvalues}{Numeric matrix (n_pcs x n_vars). Raw p-values.}
#'   \item{effect_sizes}{Numeric matrix (n_pcs x n_vars). Spearman rho for
#'     continuous variables; eta-squared for categorical variables.}
#'   \item{test_used}{Named character vector. Test applied per variable:
#'     "spearman", "kruskal", or "skipped" (degenerate variable).}
#'   \item{var_explained}{Numeric vector. Proportion of variance explained per
#'     PC (e.g. 0.25 = 25%).}
#'   \item{pca}{The \code{prcomp} object for further use or custom plotting.}
#'   \item{params}{List of parameters used.}
#' }
#'
#' @details
#' Test selection per variable:
#' \itemize{
#'   \item Continuous (> \code{categorical_threshold} unique values): Spearman
#'     rank correlation between PC scores and variable values.
#'   \item Categorical (<= \code{categorical_threshold} unique values, or
#'     character/factor): Kruskal-Wallis rank sum test across groups.
#'   \item Skipped: variables with only one unique value, or where fewer than
#'     2 complete observations per group exist.
#' }
#'
#' Effect sizes: Spearman rho (range -1 to 1); eta-squared for Kruskal-Wallis
#' (H - k + 1) / (n - k), where H is the test statistic, k is the number of
#' groups, and n is the number of observations. Negative eta-squared is set
#' to 0.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ol <- HADES_filter_olink_proteins(HADES_filter_olink(ol_raw))
#'
#' assoc <- ARTEMIS_pc_metadata_association(
#'   matrix      = ol$wide,
#'   sample_meta = ol$sample_meta,
#'   sample_col  = "SampleID",
#'   n_pcs       = 10
#' )
#'
#' print(assoc)
#' AETHER_plot_pc_association(assoc)
#' }
ARTEMIS_pc_metadata_association <- function(matrix,
                                             sample_meta,
                                             n_pcs                 = 10,
                                             ntop                  = NULL,
                                             sample_col            = NULL,
                                             metadata_cols         = NULL,
                                             scale                 = TRUE,
                                             categorical_threshold = 10L,
                                             verbose               = TRUE) {

  # ---------------------------------------------------------------------------
  # Align sample_meta to matrix columns
  # ---------------------------------------------------------------------------

  if (!is.null(sample_col)) {
    if (!sample_col %in% colnames(sample_meta)) {
      stop("sample_col '", sample_col, "' not found in sample_meta.")
    }
    rownames(sample_meta) <- as.character(sample_meta[[sample_col]])
  }

  missing_samples <- setdiff(colnames(matrix), rownames(sample_meta))
  if (length(missing_samples) > 0) {
    stop("Sample IDs in matrix missing from sample_meta: ",
         paste(head(missing_samples, 5), collapse = ", "),
         if (length(missing_samples) > 5) " ..." else "")
  }
  sample_meta <- sample_meta[colnames(matrix), , drop = FALSE]

  # ---------------------------------------------------------------------------
  # Check for NAs
  # ---------------------------------------------------------------------------

  n_na <- sum(is.na(matrix))
  if (n_na > 0) {
    warning(n_na, " NA values present in matrix. PCA requires complete data. ",
            "Filter with HADES_filter_olink_proteins() first.")
    # Proceed anyway using complete cases (proteins with no NAs)
    complete_rows <- rowSums(is.na(matrix)) == 0
    n_dropped     <- sum(!complete_rows)
    if (verbose) {
      cat("  Dropping", n_dropped, "proteins with any NA for PCA ",
          "(", sum(complete_rows), "retained)\n")
    }
    matrix <- matrix[complete_rows, , drop = FALSE]
    if (nrow(matrix) == 0) {
      stop("No complete-case proteins remain after NA removal.")
    }
  }

  # ---------------------------------------------------------------------------
  # Select metadata columns to test
  # ---------------------------------------------------------------------------

  exclude_cols <- c(sample_col)  # always exclude the ID column

  if (is.null(metadata_cols)) {
    metadata_cols <- setdiff(colnames(sample_meta), exclude_cols)
  } else {
    missing_meta <- setdiff(metadata_cols, colnames(sample_meta))
    if (length(missing_meta) > 0) {
      stop("metadata_cols not found in sample_meta: ",
           paste(missing_meta, collapse = ", "))
    }
  }

  if (length(metadata_cols) == 0) {
    stop("No metadata columns to test.")
  }

  # ---------------------------------------------------------------------------
  # Optionally subset to top N most variable features
  # ---------------------------------------------------------------------------

  if (!is.null(ntop)) {
    ntop <- as.integer(ntop)
    if (ntop >= nrow(matrix)) {
      if (verbose) cat("  ntop (", ntop, ") >= nrow(matrix) (", nrow(matrix),
                       "); using all features\n", sep = "")
    } else {
      row_vars <- apply(matrix, 1, stats::var, na.rm = TRUE)
      top_idx  <- order(row_vars, decreasing = TRUE)[seq_len(ntop)]
      matrix   <- matrix[top_idx, , drop = FALSE]
      if (verbose) cat("  Using top", ntop, "most variable features for PCA\n")
    }
  }

  # ---------------------------------------------------------------------------
  # PCA
  # ---------------------------------------------------------------------------

  if (verbose) {
    cat("Running PCA on", nrow(matrix), "features x", ncol(matrix),
        "samples...\n")
  }

  pca_result  <- stats::prcomp(t(matrix), scale. = scale, center = TRUE)
  var_exp_all <- pca_result$sdev^2 / sum(pca_result$sdev^2)

  n_pcs_actual <- min(n_pcs, ncol(pca_result$x))
  if (n_pcs_actual < n_pcs && verbose) {
    cat("  Requested", n_pcs, "PCs; only", n_pcs_actual, "available\n")
  }
  n_pcs        <- n_pcs_actual
  pc_scores    <- pca_result$x[, seq_len(n_pcs), drop = FALSE]
  var_explained <- var_exp_all[seq_len(n_pcs)]

  if (verbose) {
    cum_var <- round(100 * sum(var_explained), 1)
    cat("  PC1-", n_pcs, " explain ", cum_var, "% of variance\n", sep = "")
  }

  # ---------------------------------------------------------------------------
  # Classify metadata variables
  # ---------------------------------------------------------------------------

  test_used <- character(length(metadata_cols))
  names(test_used) <- metadata_cols

  for (col in metadata_cols) {
    vals <- sample_meta[[col]]
    n_unique <- length(unique(stats::na.omit(vals)))

    if (n_unique <= 1) {
      test_used[col] <- "skipped"
    } else if (is.character(vals) || is.factor(vals) ||
               n_unique <= categorical_threshold) {
      test_used[col] <- "kruskal"
    } else {
      test_used[col] <- "spearman"
    }
  }

  n_skipped <- sum(test_used == "skipped")
  if (verbose && n_skipped > 0) {
    cat("  Skipped", n_skipped, "degenerate variable(s):",
        paste(names(test_used)[test_used == "skipped"], collapse = ", "), "\n")
  }

  # ---------------------------------------------------------------------------
  # Run association tests
  # ---------------------------------------------------------------------------

  pc_names  <- paste0("PC", seq_len(n_pcs))
  pval_mat  <- matrix(NA_real_, nrow = n_pcs, ncol = length(metadata_cols),
                      dimnames = list(pc_names, metadata_cols))
  eff_mat   <- matrix(NA_real_, nrow = n_pcs, ncol = length(metadata_cols),
                      dimnames = list(pc_names, metadata_cols))

  for (col in metadata_cols) {
    if (test_used[col] == "skipped") next

    vals     <- sample_meta[[col]]
    complete <- !is.na(vals)

    for (pc_i in seq_len(n_pcs)) {
      scores <- pc_scores[, pc_i]

      if (test_used[col] == "spearman") {
        # Spearman correlation: PC scores vs continuous variable
        ct <- tryCatch(
          stats::cor.test(scores[complete],
                          as.numeric(vals[complete]),
                          method = "spearman",
                          exact  = FALSE),
          error = function(e) NULL
        )
        if (!is.null(ct)) {
          pval_mat[pc_i, col] <- ct$p.value
          eff_mat[pc_i,  col] <- ct$estimate  # rho
        }

      } else {
        # Kruskal-Wallis: PC scores by categorical groups
        grps <- as.character(vals[complete])
        scrs <- scores[complete]

        # Need at least 2 groups with >=1 observation each
        grp_tbl <- table(grps)
        if (sum(grp_tbl >= 1) < 2) next

        kt <- tryCatch(
          stats::kruskal.test(scrs ~ factor(grps)),
          error = function(e) NULL
        )
        if (!is.null(kt)) {
          pval_mat[pc_i, col] <- kt$p.value
          # Eta-squared: (H - k + 1) / (n - k)
          H <- kt$statistic
          k <- length(unique(grps))
          n <- length(grps)
          eta_sq <- max(0, (H - k + 1) / (n - k))
          eff_mat[pc_i, col] <- eta_sq
        }
      }
    }
  }

  if (verbose) {
    cat("  Association tests complete.\n")
    # Report top associations
    min_p    <- apply(pval_mat, 2, min, na.rm = TRUE)
    top_vars <- sort(min_p)[seq_len(min(5, length(min_p)))]
    cat("  Top associated variables (by min p across PCs):\n")
    for (nm in names(top_vars)) {
      cat("    ", nm, ": p =",
          format(top_vars[nm], digits = 2, scientific = TRUE), "\n")
    }
  }

  # ---------------------------------------------------------------------------
  # Assemble output
  # ---------------------------------------------------------------------------

  result <- list(
    pvalues      = pval_mat,
    effect_sizes = eff_mat,
    test_used    = test_used,
    var_explained = var_explained,
    pca          = pca_result,
    params       = list(
      n_pcs                 = n_pcs,
      ntop                  = ntop,
      scale                 = scale,
      categorical_threshold = categorical_threshold,
      metadata_cols         = metadata_cols
    )
  )
  class(result) <- c("artemis_pc_assoc", "list")

  return(result)
}


#' @method print artemis_pc_assoc
#' @export
print.artemis_pc_assoc <- function(x, ...) {
  cat("PC-Metadata Association\n")
  cat("-----------------------\n")
  cat("PCs tested     :", x$params$n_pcs, "\n")
  cum_var <- round(100 * sum(x$var_explained), 1)
  cat("Variance (PC1-", x$params$n_pcs, "): ", cum_var, "%\n", sep = "")

  tested   <- names(x$test_used)[x$test_used != "skipped"]
  skipped  <- names(x$test_used)[x$test_used == "skipped"]
  cat("Variables tested:", length(tested), "\n")
  if (length(skipped) > 0) {
    cat("Skipped (degenerate):", paste(skipped, collapse = ", "), "\n")
  }

  cat("\nSignificant associations (p < 0.05):\n")
  sig_any <- apply(x$pvalues, 2, function(p) any(p < 0.05, na.rm = TRUE))
  if (!any(sig_any)) {
    cat("  None\n")
  } else {
    for (col in names(sig_any)[sig_any]) {
      sig_pcs <- rownames(x$pvalues)[!is.na(x$pvalues[, col]) &
                                       x$pvalues[, col] < 0.05]
      cat("  ", col, ": ", paste(sig_pcs, collapse = ", "), "\n", sep = "")
    }
  }
  invisible(x)
}
