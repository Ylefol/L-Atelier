#' Differential Analysis for Log2-Scale Data (limma)
#'
#' @description Performs differential analysis using limma linear models on
#' any log2-scale numeric matrix. Suitable for Olink NPX data,
#' log2-transformed proteomics, or microarray expression. For count-based
#' assays (RNA-seq, ATAC-seq) use \code{ARTEMIS_differential_counts()} instead.
#'
#' @param matrix Numeric matrix, features × samples (log2-scale). Rownames
#'   are used as feature identifiers and must be present.
#' @param sample_meta Data.frame with one row per sample. Rows must correspond
#'   to columns of \code{matrix}, matched by \code{sample_col} (or rownames if
#'   \code{sample_col = NULL}).
#' @param group_col Character. Column in \code{sample_meta} defining the
#'   comparison groups.
#' @param reference Character. Reference/baseline group label (denominator in
#'   fold change). Positive logFC = higher in \code{experiment} vs \code{reference}.
#' @param experiment Character. Experimental group label (numerator in fold change).
#' @param sample_col Character or NULL. Column in \code{sample_meta} holding
#'   sample IDs matching \code{colnames(matrix)}. NULL = use
#'   \code{rownames(sample_meta)}. Default: NULL.
#' @param covariates Character vector or NULL. Column names in \code{sample_meta}
#'   to include as additive covariates in the design matrix. Default: NULL.
#' @param block_col Character or NULL. Column in \code{sample_meta} identifying
#'   the blocking variable for repeated measures (e.g. SubjectID for longitudinal
#'   data). When provided, \code{limma::duplicateCorrelation()} estimates the
#'   within-block correlation and incorporates it into the model. Default: NULL.
#' @param feature_meta Data.frame or NULL. Optional per-feature annotations to
#'   join onto the results table. Matched first by rownames (if meaningful),
#'   then by any column whose values match the matrix rownames. Default: NULL.
#' @param feature_id_col Character or NULL. Column in \code{feature_meta} to
#'   use as \code{feature_id} in the results instead of the matrix rownames.
#'   The original matrix rownames are retained as \code{original_id}.
#'   Useful for replacing internal IDs (e.g. OlinkID) with gene symbols
#'   (e.g. \code{"Assay"} from \code{ol$assay_meta}). A warning is issued if
#'   the column contains duplicates. Default: NULL.
#' @param fdr Numeric. FDR threshold used for the summary and the \code{sig}
#'   flag in results. Default: 0.05.
#' @param lfc Numeric. Minimum absolute log2 fold change for the \code{sig}
#'   flag (applied in addition to FDR). Default: 0.
#' @param verbose Logical. Print progress and summary. Default: TRUE.
#'
#' @return An S3 object of class \code{"artemis_limma"} containing:
#' \describe{
#'   \item{results}{Data.frame sorted by pvalue: feature_id, (original_id if
#'     \code{feature_id_col} used), (feature_meta columns if provided),
#'     log2FoldChange, AveExpr, t, pvalue, padj, B, sig (TRUE/FALSE).}
#'   \item{fit}{The limma MArrayLM object after \code{eBayes()} (for advanced use).}
#'   \item{summary}{List: n_features, n_tested, n_sig_up, n_sig_down, fdr, lfc.}
#'   \item{comparison}{Character string "experiment vs reference".}
#'   \item{norm_matrix}{The input matrix subset to the two groups analysed.}
#'   \item{params}{List of parameters used (for provenance).}
#' }
#'
#' @details
#' The design matrix uses a no-intercept parameterisation (\code{~ 0 + group})
#' with an explicit contrast between the two groups. Covariates are added as
#' additive terms: \code{~ 0 + group + covariate1 + ...}.
#'
#' For Olink NPX data, NPX is already log2-scale — do NOT transform the matrix
#' before passing it here. For mass spectrometry data, ensure values are
#' log2-transformed and normalised before calling this function.
#'
#' The \code{sig} column combines both FDR and LFC thresholds:
#' \code{padj < fdr & abs(log2FoldChange) >= lfc}. When \code{lfc = 0} (default),
#' only FDR filtering applies.
#'
#' Background for ORA/GSEA: all tested features are present in \code{$results}
#' (including non-significant). Use \code{results$feature_id} as the ORA
#' background, consistent with how \code{ARTEMIS_differential_counts()} is used.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Olink: two-group comparison
#' ol <- ELEUTHIA_load_olink("data.parquet", metadata_file = "layout.xlsx")
#' de <- ARTEMIS_limma_de(
#'   matrix      = ol$wide,
#'   sample_meta = ol$sample_meta,
#'   group_col   = "Group",
#'   reference   = "control",
#'   experiment  = "patient",
#'   feature_meta = ol$assay_meta
#' )
#' print(de)
#' sig_proteins <- subset(de$results, sig)
#'
#' # Longitudinal design (block by subject)
#' de_long <- ARTEMIS_limma_de(
#'   matrix      = patient_matrix,
#'   sample_meta = patient_meta,
#'   group_col   = "Timepoint",
#'   reference   = "1",
#'   experiment  = "3",
#'   block_col   = "SubjectID",
#'   feature_meta = ol$assay_meta
#' )
#'
#' # With covariate adjustment
#' de_adj <- ARTEMIS_limma_de(
#'   matrix      = ol$wide,
#'   sample_meta = ol$sample_meta,
#'   group_col   = "Group",
#'   reference   = "control",
#'   experiment  = "patient",
#'   covariates  = c("Age", "Sex")
#' )
#'
#' # Use gene symbols as feature_id instead of OlinkIDs
#' de <- ARTEMIS_limma_de(
#'   matrix        = ol$wide,
#'   sample_meta   = ol$sample_meta,
#'   group_col     = "Group",
#'   reference     = "control",
#'   experiment    = "patient",
#'   feature_meta  = ol$assay_meta,
#'   feature_id_col = "Assay"
#' )
#' }
ARTEMIS_limma_de <- function(matrix,
                              sample_meta,
                              group_col,
                              reference,
                              experiment,
                              sample_col     = NULL,
                              covariates     = NULL,
                              block_col      = NULL,
                              feature_meta   = NULL,
                              feature_id_col = NULL,
                              fdr            = 0.05,
                              lfc            = 0,
                              verbose        = TRUE) {

  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required. Install with: BiocManager::install('limma')")
  }

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
    stop("Sample IDs in matrix columns missing from sample_meta: ",
         paste(head(missing_samples, 5), collapse = ", "),
         if (length(missing_samples) > 5) " ..." else "")
  }

  # Reorder sample_meta rows to match matrix column order
  sample_meta <- sample_meta[colnames(matrix), , drop = FALSE]

  # ---------------------------------------------------------------------------
  # Validate groups
  # ---------------------------------------------------------------------------

  if (!group_col %in% colnames(sample_meta)) {
    stop("group_col '", group_col, "' not found in sample_meta.")
  }

  all_groups <- as.character(sample_meta[[group_col]])
  unique_grps <- unique(all_groups)

  if (!reference %in% unique_grps) {
    stop("reference '", reference, "' not found in group_col '", group_col,
         "'. Groups present: ", paste(unique_grps, collapse = ", "))
  }
  if (!experiment %in% unique_grps) {
    stop("experiment '", experiment, "' not found in group_col '", group_col,
         "'. Groups present: ", paste(unique_grps, collapse = ", "))
  }

  # ---------------------------------------------------------------------------
  # Subset to the two groups
  # ---------------------------------------------------------------------------

  keep_mask <- all_groups %in% c(reference, experiment)
  mat       <- matrix[, keep_mask, drop = FALSE]
  meta      <- sample_meta[keep_mask, , drop = FALSE]
  grp       <- as.character(meta[[group_col]])

  n_ref <- sum(grp == reference)
  n_exp <- sum(grp == experiment)

  if (verbose) {
    cat("[ARTEMIS] Differential analysis (limma):\n")
    cat("    Reference : ", reference,  " (n = ", n_ref, ")\n", sep = "")
    cat("    Experiment: ", experiment, " (n = ", n_exp, ")\n", sep = "")
    if (!is.null(covariates)) {
      cat("    Covariates:", paste(covariates, collapse = ", "), "\n")
    }
    if (!is.null(block_col)) {
      cat("    Block     :", block_col,
          "(repeated measures via duplicateCorrelation)\n")
    }
  }

  # ---------------------------------------------------------------------------
  # Validate optional columns
  # ---------------------------------------------------------------------------

  if (!is.null(covariates)) {
    missing_cov <- setdiff(covariates, colnames(meta))
    if (length(missing_cov) > 0) {
      stop("Covariate columns not found in sample_meta: ",
           paste(missing_cov, collapse = ", "))
    }
  }

  if (!is.null(block_col) && !block_col %in% colnames(meta)) {
    stop("block_col '", block_col, "' not found in sample_meta.")
  }

  # ---------------------------------------------------------------------------
  # Build design matrix
  # ---------------------------------------------------------------------------

  # Safe column names for contrasts (model.matrix prepends the variable name;
  # make.names ensures valid R symbols)
  safe_ref <- make.names(reference)
  safe_exp <- make.names(experiment)

  # Build a data.frame for model.matrix so formula evaluation is unambiguous
  grp_factor  <- factor(grp, levels = c(reference, experiment))
  design_data <- data.frame(grp_factor = grp_factor,
                             stringsAsFactors = FALSE)

  if (!is.null(covariates)) {
    for (cov in covariates) {
      val <- meta[[cov]]
      if (is.character(val)) val <- factor(val)
      design_data[[cov]] <- val
    }
    formula_str <- paste("~ 0 + grp_factor +", paste(covariates, collapse = " + "))
  } else {
    formula_str <- "~ 0 + grp_factor"
  }

  design     <- stats::model.matrix(stats::as.formula(formula_str),
                                    data = design_data)
  grp_cols   <- grep("^grp_factor", colnames(design))
  colnames(design)[grp_cols] <- c(safe_ref, safe_exp)

  # ---------------------------------------------------------------------------
  # Contrast
  # ---------------------------------------------------------------------------

  contrast_str <- paste0(safe_exp, " - ", safe_ref)
  contrast_mat <- limma::makeContrasts(contrasts = contrast_str, levels = design)

  # ---------------------------------------------------------------------------
  # Fit
  # ---------------------------------------------------------------------------

  if (!is.null(block_col)) {
    block_vec <- meta[[block_col]]
    if (verbose) cat("    Estimating within-block correlation...\n")
    corfit <- limma::duplicateCorrelation(mat, design, block = block_vec)
    if (verbose) {
      cat("    Consensus correlation:", round(corfit$consensus, 3), "\n")
    }
    fit <- limma::lmFit(mat, design,
                        block       = block_vec,
                        correlation = corfit$consensus)
  } else {
    fit <- limma::lmFit(mat, design)
  }

  fit2 <- limma::contrasts.fit(fit, contrast_mat)
  fit2 <- limma::eBayes(fit2)

  # ---------------------------------------------------------------------------
  # Extract results
  # ---------------------------------------------------------------------------

  tt <- limma::topTable(fit2, coef = 1, number = Inf, sort.by = "none")

  results_df <- data.frame(
    feature_id     = rownames(tt),
    log2FoldChange = tt$logFC,
    AveExpr        = tt$AveExpr,
    t              = tt$t,
    pvalue         = tt$P.Value,
    padj           = tt$adj.P.Val,
    B              = tt$B,
    stringsAsFactors = FALSE
  )

  results_df$sig <- (!is.na(results_df$padj)) &
                    (results_df$padj < fdr) &
                    (abs(results_df$log2FoldChange) >= lfc)

  # ---------------------------------------------------------------------------
  # Join feature metadata
  # ---------------------------------------------------------------------------

  if (!is.null(feature_meta) && is.data.frame(feature_meta)) {
    feature_ids <- results_df$feature_id
    fm_matched  <- NULL

    # Try 1: rownames are meaningful IDs (not just sequential integers)
    rn <- rownames(feature_meta)
    if (!is.null(rn) && any(rn %in% feature_ids) &&
        suppressWarnings(anyNA(as.integer(rn)))) {
      fm_matched <- feature_meta[feature_ids, , drop = FALSE]
      rownames(fm_matched) <- NULL

    } else {
      # Try 2: find a column whose values match the feature IDs
      for (col in colnames(feature_meta)) {
        if (any(as.character(feature_meta[[col]]) %in% feature_ids)) {
          idx        <- match(feature_ids, as.character(feature_meta[[col]]))
          fm_matched <- feature_meta[idx, , drop = FALSE]
          rownames(fm_matched) <- NULL
          break
        }
      }
    }

    if (!is.null(fm_matched)) {
      new_cols <- setdiff(colnames(fm_matched), colnames(results_df))
      if (length(new_cols) > 0) {
        results_df <- cbind(results_df, fm_matched[, new_cols, drop = FALSE])
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Apply feature_id_col: replace feature_id with a human-readable column
  # ---------------------------------------------------------------------------

  if (!is.null(feature_id_col)) {
    if (!feature_id_col %in% colnames(results_df)) {
      warning("feature_id_col '", feature_id_col, "' not found in results ",
              "after joining feature_meta — feature_id unchanged.")
    } else {
      new_ids <- as.character(results_df[[feature_id_col]])
      if (anyDuplicated(new_ids[!is.na(new_ids)])) {
        warning("feature_id_col '", feature_id_col, "' contains duplicate ",
                "values; results may be ambiguous.")
      }
      # Rename original feature_id → original_id, promote new column → feature_id
      colnames(results_df)[colnames(results_df) == "feature_id"]    <- "original_id"
      colnames(results_df)[colnames(results_df) == feature_id_col]  <- "feature_id"
      # Reorder: feature_id first, original_id second, then everything else
      other_cols <- setdiff(colnames(results_df), c("feature_id", "original_id"))
      results_df <- results_df[, c("feature_id", "original_id", other_cols)]
    }
  }

  # Sort by pvalue
  results_df <- results_df[order(results_df$pvalue), ]
  rownames(results_df) <- NULL

  # ---------------------------------------------------------------------------
  # Summary
  # ---------------------------------------------------------------------------

  n_sig_up   <- sum(results_df$sig & results_df$log2FoldChange >  0, na.rm = TRUE)
  n_sig_down <- sum(results_df$sig & results_df$log2FoldChange <  0, na.rm = TRUE)

  summary_stats <- list(
    n_features = nrow(results_df),
    n_tested   = sum(!is.na(results_df$padj)),
    n_sig_up   = n_sig_up,
    n_sig_down = n_sig_down,
    fdr        = fdr,
    lfc        = lfc
  )

  if (verbose) {
    cat("\n[ARTEMIS] Results (FDR < ", fdr, if (lfc > 0) paste0(", |log2FC| >= ", lfc) else "",
        "):\n", sep = "")
    cat("    Features tested:", summary_stats$n_tested, "\n")
    cat("    Significant    :", n_sig_up + n_sig_down,
        " (up: ", n_sig_up, ", down: ", n_sig_down, ")\n", sep = "")
    if ((n_sig_up + n_sig_down) > 0) {
      cat("    Top hits:\n")
      sig_hits <- head(results_df[results_df$sig, ], 5)
      for (i in seq_len(nrow(sig_hits))) {
        dir <- if (sig_hits$log2FoldChange[i] > 0) "UP" else "DOWN"
        cat("        ", sig_hits$feature_id[i],
            ": log2FC = ", round(sig_hits$log2FoldChange[i], 2),
            " (", dir, "), padj = ",
            format(sig_hits$padj[i], digits = 2, scientific = TRUE),
            "\n", sep = "")
      }
    }
    cat("\n")
  }

  # ---------------------------------------------------------------------------
  # Assemble output
  # ---------------------------------------------------------------------------

  result <- list(
    results     = results_df,
    fit         = fit2,
    summary     = summary_stats,
    comparison  = paste0(experiment, " vs ", reference),
    norm_matrix = mat,
    params      = list(
      group_col      = group_col,
      reference      = reference,
      experiment     = experiment,
      covariates     = covariates,
      block_col      = block_col,
      feature_id_col = feature_id_col,
      fdr            = fdr,
      lfc            = lfc
    )
  )
  class(result) <- c("artemis_limma", "list")

  return(result)
}


#' @method print artemis_limma
#' @export
print.artemis_limma <- function(x, ...) {
  cat("limma DE Result\n")
  cat("------------------------------\n")
  cat("Comparison : ", x$comparison, "\n", sep = "")
  cat("Features   : ", x$summary$n_features, "\n", sep = "")
  n_sig <- x$summary$n_sig_up + x$summary$n_sig_down
  cat("Significant (FDR < ", x$params$fdr,
      if (x$params$lfc > 0) paste0(", |logFC| >= ", x$params$lfc) else "",
      "): ", n_sig,
      " (up: ", x$summary$n_sig_up,
      ", down: ", x$summary$n_sig_down, ")\n", sep = "")
  if (!is.null(x$params$covariates)) {
    cat("Covariates : ", paste(x$params$covariates, collapse = ", "), "\n", sep = "")
  }
  if (!is.null(x$params$block_col)) {
    cat("Block: ", x$params$block_col, "\n", sep = "")
  }
  cat("\nSlots: $results, $fit, $summary, $comparison, $norm_matrix, $params\n")
  invisible(x)
}
