# GAIA/Demeter/correlation_utils.R
# Correlation and p-value utilities
#
# General-purpose correlation functions usable across modules.


#' Calculate correlation matrix with p-values
#'
#' Computes pairwise correlations between all columns of a data.frame/matrix,
#' along with p-values for each correlation.
#'
#' @param data Data.frame or matrix with variables as columns.
#' @param method Correlation method: "spearman", "pearson", or "kendall".
#'   Default: "spearman".
#' @param alternative Alternative hypothesis: "two.sided", "less", "greater".
#'   Default: "two.sided".
#' @param adjust P-value adjustment method: "BH", "bonferroni", "none", etc.
#'   Default: "none".
#' @param use How to handle missing values: "pairwise.complete.obs" (default),
#'   "complete.obs", "everything".
#'
#' @return A list with:
#'   \item{correlation}{Correlation matrix}
#'   \item{pvalue}{P-value matrix}
#'   \item{padj}{Adjusted p-value matrix (if adjust != "none")}
#'   \item{n}{Sample size matrix (pairwise)}
#'   \item{method}{Correlation method used}
#'
#' @examples
#' \dontrun{
#' data <- data.frame(a = rnorm(100), b = rnorm(100), c = rnorm(100))
#' result <- DEMETER_correlation_matrix(data)
#' result$correlation
#' result$pvalue
#'
#' }
#' @export
DEMETER_correlation_matrix <- function(data,
                                        method = "spearman",
                                        alternative = "two.sided",
                                        adjust = "none",
                                        use = "pairwise.complete.obs") {

  data <- as.matrix(data)
  n_vars <- ncol(data)
  var_names <- colnames(data)

  if (is.null(var_names)) {
    var_names <- paste0("V", seq_len(n_vars))
  }

  # Initialize matrices
  cor_mat <- matrix(NA, n_vars, n_vars)
  p_mat <- matrix(NA, n_vars, n_vars)
  n_mat <- matrix(NA, n_vars, n_vars)

  rownames(cor_mat) <- colnames(cor_mat) <- var_names
  rownames(p_mat) <- colnames(p_mat) <- var_names
  rownames(n_mat) <- colnames(n_mat) <- var_names

  # Fill diagonal
  diag(cor_mat) <- 1
  diag(p_mat) <- 0

  # Calculate pairwise correlations
  for (i in seq_len(n_vars - 1)) {
    for (j in (i + 1):n_vars) {
      x <- data[, i]
      y <- data[, j]

      # Remove NAs pairwise
      complete <- complete.cases(x, y)
      n_mat[i, j] <- n_mat[j, i] <- sum(complete)

      if (sum(complete) >= 3) {
        test <- cor.test(x[complete], y[complete],
                         method = method,
                         alternative = alternative,
                         exact = FALSE)

        cor_mat[i, j] <- cor_mat[j, i] <- test$estimate
        p_mat[i, j] <- p_mat[j, i] <- test$p.value
      }
    }
  }

  diag(n_mat) <- colSums(!is.na(data))

  # Adjust p-values if requested
  if (adjust != "none") {
    # Get upper triangle p-values (avoid double-counting)
    upper_p <- p_mat[upper.tri(p_mat)]
    upper_padj <- p.adjust(upper_p, method = adjust)

    padj_mat <- matrix(NA, n_vars, n_vars)
    rownames(padj_mat) <- colnames(padj_mat) <- var_names
    diag(padj_mat) <- 0
    padj_mat[upper.tri(padj_mat)] <- upper_padj
    padj_mat[lower.tri(padj_mat)] <- t(padj_mat)[lower.tri(padj_mat)]
  } else {
    padj_mat <- NULL
  }

  result <- list(
    correlation = cor_mat,
    pvalue = p_mat,
    padj = padj_mat,
    n = n_mat,
    method = method,
    alternative = alternative,
    adjust = adjust
  )

  return(result)
}


#' Format p-values for display
#'
#' Converts p-values to formatted strings, optionally in scientific notation.
#' Non-significant values can be replaced with empty strings or NA.
#'
#' @param pvalues Numeric vector or matrix of p-values.
#' @param format Output format: "scientific" (e.g., "1.23e-04"), "stars"
#'   (*, **, ***), "decimal" (e.g., "0.001"), or "hybrid" (decimal if >= 0.001,
#'   scientific otherwise). Default: "scientific".
#' @param digits Number of significant digits. Default: 2.
#' @param sig_threshold P-value threshold for significance. Values above this
#'   become empty strings (if hide_ns = TRUE). Default: 0.05.
#' @param hide_ns Logical. Replace non-significant p-values with empty string.
#'   Default: FALSE.
#' @param ns_string String to use for non-significant values when hide_ns = TRUE.
#'   Default: "".
#'
#' @return Character vector or matrix of formatted p-values (same shape as input).
#'
#' @examples
#' \dontrun{
#' pvals <- c(0.001, 0.023, 0.15, 0.0001)
#' DEMETER_format_pvalues(pvals)
#' DEMETER_format_pvalues(pvals, format = "stars")
#' DEMETER_format_pvalues(pvals, hide_ns = TRUE)
#'
#' }
#' @export
DEMETER_format_pvalues <- function(pvalues,
                                    format = "scientific",
                                    digits = 2,
                                    sig_threshold = 0.05,
                                    hide_ns = FALSE,
                                    ns_string = "") {

  # Store original dimensions
  orig_dim <- dim(pvalues)
  orig_dimnames <- dimnames(pvalues)
  pvalues <- as.vector(pvalues)

  formatted <- character(length(pvalues))

  for (i in seq_along(pvalues)) {
    p <- pvalues[i]

    if (is.na(p)) {
      formatted[i] <- ""
      next
    }

    # Check significance
    if (hide_ns && p >= sig_threshold) {
      formatted[i] <- ns_string
      next
    }

    # Format based on type
    formatted[i] <- switch(
      format,
      "scientific" = formatC(p, format = "e", digits = digits),
      "decimal" = formatC(signif(p, digits + 1), format = "fg", digits = digits),
      "stars" = {
        if (p < 0.001) "***"
        else if (p < 0.01) "**"
        else if (p < 0.05) "*"
        else if (p < 0.1) "."
        else ""
      },
      "hybrid" = {
        if (p >= 0.001) formatC(signif(p, digits + 1), format = "fg", digits = digits)
        else formatC(p, format = "e", digits = digits)
      },
      formatC(p, format = "e", digits = digits)  # Default to scientific
    )
  }

  # Restore dimensions if matrix
  if (!is.null(orig_dim)) {
    formatted <- matrix(formatted, nrow = orig_dim[1], ncol = orig_dim[2])
    dimnames(formatted) <- orig_dimnames
  }

  return(formatted)
}


#' Create lower/upper triangular p-value matrix for plotting
#'
#' Formats a p-value matrix for use with correlation plots, where only
#' the upper or lower triangle is shown.
#'
#' @param pvalue_matrix Square matrix of p-values.
#' @param triangle Which triangle to keep: "upper" or "lower". Default: "upper".
#' @param format Format for p-values (passed to DEMETER_format_pvalues).
#' @param sig_threshold Significance threshold.
#' @param hide_ns Hide non-significant values.
#'
#' @return Character matrix with formatted p-values in specified triangle,
#'   empty strings elsewhere.
#'
#' @export
DEMETER_format_pvalue_triangle <- function(pvalue_matrix,
                                            triangle = "upper",
                                            format = "scientific",
                                            sig_threshold = 0.05,
                                            hide_ns = TRUE) {

  formatted <- DEMETER_format_pvalues(
    pvalue_matrix,
    format = format,
    sig_threshold = sig_threshold,
    hide_ns = hide_ns
  )

  # Blank out the opposite triangle and diagonal
  if (triangle == "upper") {
    formatted[lower.tri(formatted, diag = TRUE)] <- ""
  } else {
    formatted[upper.tri(formatted, diag = TRUE)] <- ""
  }

  return(formatted)
}
