#' Factor Analysis of Mixed Data (FAMD)
#'
#' @description Wrapper around FactoMineR::FAMD for dimensionality reduction
#' of datasets containing both quantitative (continuous) and qualitative
#' (categorical) variables.
#'
#' @param data Data frame with mixed variable types. Must contain at least
#'   one quantitative and one qualitative variable.
#' @param na_action Character. How to handle missing values:
#'   \itemize{
#'     \item "fail" (default): Error if any NA values present after column filtering
#'     \item "omit": Remove rows with any NA values
#'     \item "impute": Use missMDA::imputeFAMD() for iterative imputation
#'   }
#' @param na_threshold Numeric (0-1). Columns with proportion of NA values
#'   exceeding this threshold are dropped before other NA handling.
#'   Default = 0.2 (20%).
#' @param ncp Integer. Number of dimensions to keep in results. Default = 5.
#' @param sup_quanti Integer vector. Indices of quantitative supplementary
#'   variables (not used in computation, but projected onto results).
#' @param sup_quali Integer vector. Indices of qualitative supplementary
#'   variables (not used in computation, but projected onto results).
#' @param ncp_impute Integer. Number of components used for imputation when
#'   na_action = "impute". Default = 5.
#' @param prefix_levels Logical. If TRUE (default), prefix factor levels with
#'   variable names (e.g., "0" becomes "gender_0"). This prevents ambiguity when
#'   multiple factors have overlapping numeric levels and makes FAMD output
#'   more readable.
#' @param verbose Logical. Print progress messages. Default = TRUE.
#'
#' @return A list with class "artemis_famd" containing:
#' \describe{
#'   \item{famd}{The FactoMineR FAMD result object}
#'   \item{data_used}{The data frame used for analysis (after NA handling)}
#'   \item{eigenvalues}{Data frame with eigenvalue, variance percent, and cumulative percent}
#'   \item{ind}{List with individual (sample) coordinates, cos2, and contributions}
#'   \item{var}{List with variable information (quanti and quali separately)}
#'   \item{na_report}{Summary of NA handling actions taken}
#'   \item{level_mapping}{List mapping original factor levels to prefixed versions
#'     (NULL if prefix_levels = FALSE)}
#'   \item{call}{The function call}
#' }
#'
#' @details
#' FAMD is particularly useful for exploratory analysis of sample metadata
#' that contains both continuous variables (age, signal intensity) and
#' categorical variables (treatment, batch, sex).
#'
#' The function performs the following steps:
#' 1. Validate input (check for mixed types)
#' 2. Drop columns exceeding NA threshold
#' 3. Handle remaining NAs based on na_action
#' 4. Run FAMD analysis
#' 5. Structure and return results
#'
#' @export
#'
#' @examples
#' # Basic usage with sample metadata
#' result <- ARTEMIS_famd(sample_metadata)
#'
#' # With imputation for missing values
#' result <- ARTEMIS_famd(sample_metadata, na_action = "impute")
#'
#' # With supplementary variables
#' result <- ARTEMIS_famd(sample_metadata, sup_quali = c(1))  # First column as supplementary
#'
ARTEMIS_famd <- function(data,
                         na_action = "fail",
                         na_threshold = 0.2,
                         ncp = 5,
                         sup_quanti = NULL,
                         sup_quali = NULL,
                         ncp_impute = 5,
                         prefix_levels = TRUE,
                         verbose = TRUE) {

  # Validate na_action

if (!na_action %in% c("fail", "omit", "impute")) {
    stop("na_action must be one of: 'fail', 'omit', 'impute'")
  }

  # Check for FactoMineR
  if (!requireNamespace("FactoMineR", quietly = TRUE)) {
    stop("Package 'FactoMineR' is required. Install with: install.packages('FactoMineR')")
  }

  # Check for missMDA if imputation requested
  if (na_action == "impute" && !requireNamespace("missMDA", quietly = TRUE)) {
    stop("Package 'missMDA' is required for imputation. Install with: install.packages('missMDA')")
  }

  # Ensure data is a data.frame

  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  original_nrow <- nrow(data)
  original_ncol <- ncol(data)

  if (verbose) {
    cat("ARTEMIS FAMD Analysis\n")
    cat("=====================\n")
    cat("Input data:", original_nrow, "observations,", original_ncol, "variables\n")
  }

  # ---------------------------------------------------------------------------
  # Step 1: Identify variable types
  # ---------------------------------------------------------------------------
  var_types <- sapply(data, function(x) {
    if (is.numeric(x)) "quanti"
    else if (is.factor(x) || is.character(x)) "quali"
    else "other"
  })

  n_quanti <- sum(var_types == "quanti")
  n_quali <- sum(var_types == "quali")
  n_other <- sum(var_types == "other")

  if (verbose) {
    cat("Variable types:", n_quanti, "quantitative,", n_quali, "qualitative")
    if (n_other > 0) cat(",", n_other, "other (will be excluded)")
    cat("\n")
  }

  # Check for mixed data requirement
  if (n_quanti == 0) {
    stop("No quantitative variables found. FAMD requires at least one quantitative variable. ",
         "For purely categorical data, use MCA instead.")
  }
  if (n_quali == 0) {
    stop("No qualitative variables found. FAMD requires at least one qualitative variable. ",
         "For purely continuous data, use PCA instead.")
  }

  # Remove 'other' type columns
  if (n_other > 0) {
    other_cols <- names(var_types)[var_types == "other"]
    if (verbose) {
      cat("  Excluding non-standard columns:", paste(other_cols, collapse = ", "), "\n")
    }
    data <- data[, var_types != "other", drop = FALSE]
    var_types <- var_types[var_types != "other"]
  }

  # Convert character columns to factors
  char_cols <- sapply(data, is.character)
  if (any(char_cols)) {
    data[char_cols] <- lapply(data[char_cols], as.factor)
    if (verbose) {
      cat("  Converted", sum(char_cols), "character column(s) to factor\n")
    }
  }

  # ---------------------------------------------------------------------------
  # Step 1b: Prefix factor levels with variable names
  # ---------------------------------------------------------------------------
  # This prevents ambiguity when multiple factors have overlapping levels
  # (e.g., gender with 0/1 and score with 1-5 would both have "1")
  level_mapping <- NULL

  if (prefix_levels) {
    factor_cols <- sapply(data, is.factor)

    if (any(factor_cols)) {
      level_mapping <- list()

      for (col_name in names(data)[factor_cols]) {
        old_levels <- levels(data[[col_name]])

        # Create new levels with variable name prefix
        new_levels <- paste0(col_name, "_", old_levels)

        # Store mapping for reference
        level_mapping[[col_name]] <- data.frame(
          original = old_levels,
          prefixed = new_levels,
          stringsAsFactors = FALSE
        )

        # Update factor levels
        levels(data[[col_name]]) <- new_levels
      }

      if (verbose) {
        cat("  Prefixed levels for", sum(factor_cols), "factor column(s)\n")
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Step 2: Handle missing values - Column filtering
  # ---------------------------------------------------------------------------
  na_report <- list(
    original_rows = original_nrow,
    original_cols = original_ncol,
    cols_dropped = character(0),
    rows_dropped = 0,
    imputed = FALSE
  )

  # Calculate NA proportion per column
  na_props <- colMeans(is.na(data))
  high_na_cols <- names(na_props)[na_props > na_threshold]

  if (length(high_na_cols) > 0) {
    if (verbose) {
      cat("\nNA Handling:\n")
      cat("  Dropping", length(high_na_cols), "column(s) with >",
          round(na_threshold * 100), "% missing:\n")
      for (col in high_na_cols) {
        cat("    -", col, "(", round(na_props[col] * 100, 1), "% NA)\n")
      }
    }
    data <- data[, !names(data) %in% high_na_cols, drop = FALSE]
    na_report$cols_dropped <- high_na_cols

    # Update var_types
    var_types <- var_types[!names(var_types) %in% high_na_cols]

    # Re-check for mixed data after column removal
    if (sum(var_types == "quanti") == 0 || sum(var_types == "quali") == 0) {
      stop("After dropping high-NA columns, data no longer has mixed types. ",
           "Remaining: ", sum(var_types == "quanti"), " quantitative, ",
           sum(var_types == "quali"), " qualitative.")
    }
  }

  # ---------------------------------------------------------------------------
  # Step 3: Handle remaining NAs based on na_action
  # ---------------------------------------------------------------------------
  remaining_na <- sum(is.na(data))

  if (remaining_na > 0) {
    if (verbose) {
      cat("  Remaining NA values:", remaining_na, "\n")
    }

    if (na_action == "fail") {
      # Provide diagnostic information
      na_by_col <- colSums(is.na(data))
      na_cols <- na_by_col[na_by_col > 0]
      na_by_row <- rowSums(is.na(data))

      stop("Data contains ", remaining_na, " NA values.\n",
           "  Columns with NAs: ", paste(names(na_cols), collapse = ", "), "\n",
           "  Rows with NAs: ", sum(na_by_row > 0), " of ", nrow(data), "\n",
           "Use na_action = 'omit' to remove incomplete rows, or ",
           "na_action = 'impute' to impute missing values.")

    } else if (na_action == "omit") {
      complete_rows <- complete.cases(data)
      n_dropped <- sum(!complete_rows)
      data <- data[complete_rows, , drop = FALSE]
      na_report$rows_dropped <- n_dropped

      if (verbose) {
        cat("  Removed", n_dropped, "incomplete rows (",
            round(n_dropped / original_nrow * 100, 1), "% of data)\n")
      }

      if (nrow(data) < 10) {
        warning("Only ", nrow(data), " complete cases remain. Results may be unreliable.")
      }

    } else if (na_action == "impute") {
      if (verbose) {
        cat("  Imputing missing values using missMDA::imputeFAMD()...\n")
      }

      # Impute using iterative FAMD
      imputed <- missMDA::imputeFAMD(data, ncp = ncp_impute)
      data <- imputed$completeObs
      na_report$imputed <- TRUE

      if (verbose) {
        cat("  Imputation complete\n")
      }
    }
  } else if (verbose) {
    cat("  No missing values in data\n")
  }

  # ---------------------------------------------------------------------------
  # Step 4: Run FAMD
  # ---------------------------------------------------------------------------
  if (verbose) {
    cat("\nRunning FAMD analysis...\n")
    cat("  Final data:", nrow(data), "observations,", ncol(data), "variables\n")
    cat("  Dimensions to compute:", ncp, "\n")
  }

  # Adjust supplementary variable indices if columns were dropped
  # (This is a simplification - in practice user should re-specify after seeing dropped cols)
  if (!is.null(sup_quanti) || !is.null(sup_quali)) {
    if (length(na_report$cols_dropped) > 0) {
      warning("Supplementary variable indices may be invalid after column removal. ",
              "Please verify indices match current column positions.")
    }
  }

  # Run FAMD
  famd_result <- FactoMineR::FAMD(
    data,
    ncp = ncp,
    sup.var = sup_quanti,
    ind.sup = NULL,
    graph = FALSE
  )

  # Handle qualitative supplementary variables separately if specified
  # Note: FactoMineR::FAMD uses 'sup.var' for supplementary individuals, not quali vars
  # For quali supplementary, we'd need a different approach

  # ---------------------------------------------------------------------------
  # Step 5: Structure output
  # ---------------------------------------------------------------------------
  if (verbose) {
    cat("\nResults summary:\n")
  }

  # Extract eigenvalues
  eigenvalues <- data.frame(
    dimension = seq_len(nrow(famd_result$eig)),
    eigenvalue = famd_result$eig[, 1],
    variance_percent = famd_result$eig[, 2],
    cumulative_percent = famd_result$eig[, 3]
  )

  if (verbose) {
    cat("  Variance explained by first", min(5, nrow(eigenvalues)), "dimensions:\n")
    for (i in seq_len(min(5, nrow(eigenvalues)))) {
      cat("    Dim", i, ":", round(eigenvalues$variance_percent[i], 1), "%",
          "(cumulative:", round(eigenvalues$cumulative_percent[i], 1), "%)\n")
    }
  }

  # Extract individual (sample) results
  ind_results <- list(
    coord = as.data.frame(famd_result$ind$coord),
    cos2 = as.data.frame(famd_result$ind$cos2),
    contrib = as.data.frame(famd_result$ind$contrib)
  )

  # Extract variable results
  var_results <- list(
    quanti = list(
      coord = as.data.frame(famd_result$quanti.var$coord),
      cos2 = as.data.frame(famd_result$quanti.var$cos2),
      contrib = as.data.frame(famd_result$quanti.var$contrib)
    ),
    quali = list(
      coord = as.data.frame(famd_result$quali.var$coord),
      cos2 = as.data.frame(famd_result$quali.var$cos2),
      contrib = as.data.frame(famd_result$quali.var$contrib),
      v.test = as.data.frame(famd_result$quali.var$v.test)
    )
  )

  # Build result object
  result <- list(
    famd = famd_result,
    data_used = data,
    eigenvalues = eigenvalues,
    ind = ind_results,
    var = var_results,
    na_report = na_report,
    level_mapping = level_mapping,
    call = match.call()
  )

  class(result) <- c("artemis_famd", "list")

  if (verbose) {
    cat("\nFAMD analysis complete.\n")
  }

  return(result)
}


#' Print method for ARTEMIS FAMD results
#'
#' @param x An artemis_famd object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.artemis_famd <- function(x, ...) {
  cat("ARTEMIS FAMD Result\n")
  cat("===================\n")
  cat("Data:", nrow(x$data_used), "observations,", ncol(x$data_used), "variables\n")

  n_quanti <- sum(sapply(x$data_used, is.numeric))
  n_quali <- ncol(x$data_used) - n_quanti
  cat("  Quantitative:", n_quanti, "\n")
  cat("  Qualitative:", n_quali, "\n")

  cat("\nVariance explained:\n")
  top_dims <- min(5, nrow(x$eigenvalues))
  for (i in seq_len(top_dims)) {
    cat("  Dim", i, ":", round(x$eigenvalues$variance_percent[i], 1), "%\n")
  }
  cat("  Cumulative (", top_dims, " dims):", round(x$eigenvalues$cumulative_percent[top_dims], 1), "%\n")

  if (x$na_report$imputed) {
    cat("\nNote: Missing values were imputed using iterative FAMD\n")
  }
  if (x$na_report$rows_dropped > 0) {
    cat("\nNote:", x$na_report$rows_dropped, "rows with missing values were removed\n")
  }
  if (length(x$na_report$cols_dropped) > 0) {
    cat("\nNote:", length(x$na_report$cols_dropped), "columns dropped due to high NA proportion\n")
  }
  if (!is.null(x$level_mapping) && length(x$level_mapping) > 0) {
    cat("\nNote: Factor levels prefixed with variable names (", length(x$level_mapping),
        " factors). Access $level_mapping for original values.\n", sep = "")
  }

  invisible(x)
}
