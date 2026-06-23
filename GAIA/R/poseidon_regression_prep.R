###############################################################################
###################### Regression / Group LASSO Data Prep ####################
###############################################################################
# Encoding and scaling of mixed (categorical + quantitative) data ahead of
# penalized regression (Group LASSO via grpreg) in MINERVA.

#' Encode Mixed Data for Regression Analysis
#'
#' @description Prepares a mixed (categorical + quantitative) dataset for
#' regression analysis by separating the target variable, encoding categorical
#' predictors, and creating group mappings for Group LASSO.
#'
#' @param data A data.frame with mixed data types.
#' @param target_var Character string. Name of the target variable column.
#' @param encoding Character string. Encoding method for categorical variables:
#'   \itemize{
#'     \item "dummy" (default): Treatment/dummy coding. Creates k-1 binary
#'       columns for k levels, with one level as reference.
#'     \item "effect": Effect coding. Similar to dummy but reference level
#'       coded as -1 instead of 0 (compares to grand mean).
#'   }
#' @param reference Character string or named list. How to choose reference level:
#'   \itemize{
#'     \item "first" (default): First level alphabetically
#'     \item "last": Last level alphabetically
#'     \item Named list: Specific reference per variable,
#'       e.g., list(sex = "male", treatment = "placebo")
#'   }
#' @param ordinal_vars Character vector. Names of variables that are ordinal
#'   (ordered categorical). These can be treated differently from nominal
#'   categoricals.
#' @param ordinal_treatment Character string. How to treat ordinal variables:
#'   \itemize{
#'     \item "categorical" (default): Treat as nominal (dummy coding)
#'     \item "continuous": Treat as numeric (1, 2, 3, ...)
#'     \item "polynomial": Use polynomial contrasts (linear, quadratic, etc.)
#'   }
#' @param drop_target_na Logical. Remove rows where target is NA (default = TRUE).
#' @param verbose Logical. Print encoding summary (default = TRUE).
#'
#' @return A list containing:
#' \describe{
#'   \item{X}{Numeric matrix of encoded predictors (samples x features)}
#'   \item{y}{Vector of target variable values}
#'   \item{var_mapping}{List with encoding information:
#'     \itemize{
#'       \item original_vars: Original variable names
#'       \item var_types: Type of each original variable (quantitative/categorical/ordinal)
#'       \item groups: Integer vector mapping each column in X to its variable group
#'         (for Group LASSO - all columns from same variable share group number)
#'       \item group_names: Variable name for each group number
#'       \item col_to_var: Maps each X column name to original variable
#'       \item reference_levels: Reference level for each categorical variable
#'       \item encoding: Encoding method used
#'     }
#'   }
#'   \item{sample_ids}{Row identifiers from original data}
#'   \item{target_name}{Name of target variable}
#'   \item{n_removed_na}{Number of rows removed due to NA in target}
#' }
#'
#' @details
#' This function is the first step in a supervised variable selection workflow.
#' The output is designed to work with Group LASSO (via grpreg package), where
#' the `groups` vector ensures that all dummy columns from a single categorical
#' variable are treated as a group (either all in or all out).
#'
#' For ordinal variables (like severity scores 1-9), the choice of treatment
#' affects interpretation:
#' \itemize{
#'   \item "categorical": Most flexible, allows non-linear effects, each level
#'     can have independent effect (recommended for selection stage)
#'   \item "continuous": Assumes linear effect across ordered levels
#'   \item "polynomial": Captures non-linear ordered effects with fewer parameters
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic encoding
#' encoded <- POSEIDON_encode_for_regression(
#'   data = clinical_data,
#'   target_var = "outcome"
#' )
#'
#' # With specific reference levels
#' encoded <- POSEIDON_encode_for_regression(
#'   data = clinical_data,
#'   target_var = "response",
#'   reference = list(treatment = "placebo", sex = "female")
#' )
#'
#' # With ordinal variables treated as continuous
#' encoded <- POSEIDON_encode_for_regression(
#'   data = clinical_data,
#'   target_var = "survival",
#'   ordinal_vars = c("severity_score", "stage"),
#'   ordinal_treatment = "continuous"
#' )
#'
#' }
POSEIDON_encode_for_regression <- function(data,
                                            target_var,
                                            encoding = "dummy",
                                            reference = "first",
                                            ordinal_vars = NULL,
                                            ordinal_treatment = "categorical",
                                            drop_target_na = TRUE,
                                            verbose = TRUE) {

  # Validate inputs
  if (!is.data.frame(data)) {
    stop("data must be a data.frame")
  }

  if (!target_var %in% colnames(data)) {
    stop("target_var '", target_var, "' not found in data. ",
         "Available columns: ", paste(head(colnames(data), 10), collapse = ", "),
         if (ncol(data) > 10) "..." else "")
  }

  encoding <- tolower(encoding)
  if (!encoding %in% c("dummy", "effect")) {
    stop("encoding must be 'dummy' or 'effect'")
  }

  ordinal_treatment <- tolower(ordinal_treatment)
  if (!ordinal_treatment %in% c("categorical", "continuous", "polynomial")) {
    stop("ordinal_treatment must be 'categorical', 'continuous', or 'polynomial'")
  }

  # Store sample IDs
  if (!is.null(rownames(data))) {
    sample_ids <- rownames(data)
  } else {
    sample_ids <- as.character(seq_len(nrow(data)))
  }

  # Separate target and predictors
  y <- data[[target_var]]
  predictor_data <- data[, setdiff(colnames(data), target_var), drop = FALSE]

  # Handle NA in target
  n_removed_na <- 0
  if (drop_target_na && any(is.na(y))) {
    na_idx <- is.na(y)
    n_removed_na <- sum(na_idx)
    y <- y[!na_idx]
    predictor_data <- predictor_data[!na_idx, , drop = FALSE]
    sample_ids <- sample_ids[!na_idx]

    if (verbose) {
      cat("[POSEIDON] Removed", n_removed_na, "rows with NA in target variable\n")
    }
  }

  # Classify variables
  var_names <- colnames(predictor_data)
  var_types <- character(length(var_names))
  names(var_types) <- var_names

  for (var in var_names) {
    col <- predictor_data[[var]]
    if (var %in% ordinal_vars) {
      var_types[var] <- "ordinal"
    } else if (is.factor(col) || is.character(col) || is.logical(col)) {
      var_types[var] <- "categorical"
    } else if (is.numeric(col)) {
      var_types[var] <- "quantitative"
    } else {
      var_types[var] <- "categorical"  # Default to categorical for unknown types
      if (verbose) {
        warning("Variable '", var, "' has unknown type, treating as categorical")
      }
    }
  }

  if (verbose) {
    cat("\n[POSEIDON] Variable classification:\n")
    cat("  Quantitative:", sum(var_types == "quantitative"), "\n")
    cat("  Categorical:", sum(var_types == "categorical"), "\n")
    cat("  Ordinal:", sum(var_types == "ordinal"), "\n")
  }

  # Process reference levels
  if (is.list(reference)) {
    ref_levels <- reference
  } else {
    ref_levels <- list()
  }
  ref_rule <- if (is.character(reference) && length(reference) == 1) reference else "first"

  # Build encoded matrix column by column
  encoded_cols <- list()
  groups <- integer(0)
  group_names <- character(0)
  col_to_var <- character(0)
  reference_levels <- list()
  group_counter <- 0

  for (var in var_names) {
    col <- predictor_data[[var]]
    var_type <- var_types[var]
    group_counter <- group_counter + 1

    if (var_type == "quantitative") {
      # Numeric variable: use as-is
      encoded_cols[[var]] <- as.numeric(col)
      names(encoded_cols)[length(encoded_cols)] <- var
      groups <- c(groups, group_counter)
      col_to_var <- c(col_to_var, var)

    } else if (var_type == "ordinal" && ordinal_treatment == "continuous") {
      # Ordinal treated as continuous: convert to numeric ranks
      if (is.factor(col)) {
        encoded_cols[[var]] <- as.numeric(col)
      } else {
        col_factor <- factor(col)
        encoded_cols[[var]] <- as.numeric(col_factor)
      }
      names(encoded_cols)[length(encoded_cols)] <- var
      groups <- c(groups, group_counter)
      col_to_var <- c(col_to_var, var)

    } else if (var_type == "ordinal" && ordinal_treatment == "polynomial") {
      # Ordinal with polynomial contrasts
      col_factor <- factor(col)
      n_levels <- nlevels(col_factor)

      if (n_levels < 2) {
        if (verbose) warning("Variable '", var, "' has fewer than 2 levels, skipping")
        group_counter <- group_counter - 1
        next
      }

      contrasts(col_factor) <- contr.poly(n_levels)
      mm <- model.matrix(~ col_factor - 1)
      # Remove intercept-equivalent column
      mm <- mm[, -1, drop = FALSE]

      for (j in seq_len(ncol(mm))) {
        col_name <- paste0(var, "_poly", j)
        encoded_cols[[col_name]] <- mm[, j]
        groups <- c(groups, group_counter)
        col_to_var <- c(col_to_var, var)
      }

    } else {
      # Categorical (or ordinal treated as categorical): dummy coding
      col_factor <- factor(col)
      levels_var <- levels(col_factor)
      n_levels <- length(levels_var)

      if (n_levels < 2) {
        if (verbose) warning("Variable '", var, "' has fewer than 2 levels, skipping")
        group_counter <- group_counter - 1
        next
      }

      # Determine reference level
      if (var %in% names(ref_levels)) {
        ref <- ref_levels[[var]]
        if (!ref %in% levels_var) {
          warning("Reference level '", ref, "' not found for variable '", var,
                  "'. Using first level instead.")
          ref <- levels_var[1]
        }
      } else if (ref_rule == "last") {
        ref <- levels_var[n_levels]
      } else {
        ref <- levels_var[1]
      }

      reference_levels[[var]] <- ref

      # Reorder so reference is first (will be dropped in dummy coding)
      levels_reordered <- c(ref, setdiff(levels_var, ref))
      col_factor <- factor(col_factor, levels = levels_reordered)

      # Create dummy columns
      for (lev in levels_reordered[-1]) {  # Skip reference level
        col_name <- paste0(var, "_", lev)
        if (encoding == "dummy") {
          encoded_cols[[col_name]] <- as.numeric(col_factor == lev)
        } else {  # effect coding
          dummy <- as.numeric(col_factor == lev)
          dummy[col_factor == ref] <- -1
          encoded_cols[[col_name]] <- dummy
        }
        groups <- c(groups, group_counter)
        col_to_var <- c(col_to_var, var)
      }
    }

    group_names <- c(group_names, var)
  }

  # Combine into matrix
  X <- do.call(cbind, encoded_cols)
  rownames(X) <- sample_ids
  names(groups) <- colnames(X)

  if (verbose) {
    cat("\n[POSEIDON] Encoding summary:\n")
    cat("  Original variables:", length(var_names), "\n")
    cat("  Encoded columns:", ncol(X), "\n")
    cat("  Samples:", nrow(X), "\n")
    cat("  Groups for Group LASSO:", max(groups), "\n")

    cat("\n[POSEIDON] Column breakdown by variable:\n")
    for (g in seq_along(group_names)) {
      n_cols <- sum(groups == g)
      cat("  ", group_names[g], ": ", n_cols, " column(s)\n", sep = "")
    }
  }

  # Build var_mapping
  var_mapping <- list(
    original_vars = var_names,
    var_types = var_types,
    groups = groups,
    group_names = group_names,
    col_to_var = col_to_var,
    reference_levels = reference_levels,
    encoding = encoding,
    ordinal_treatment = ordinal_treatment
  )

  return(list(
    X = X,
    y = y,
    var_mapping = var_mapping,
    sample_ids = sample_ids,
    target_name = target_var,
    n_removed_na = n_removed_na
  ))
}


#' Scale Predictors for Penalized Regression
#'
#' @description Scales predictor matrix for use in penalized regression methods.
#' Penalized regression (LASSO, Group LASSO, elastic net) requires scaling
#' because penalties are applied uniformly - without scaling, variables on
#' larger scales would be penalized less.
#'
#' @param X Numeric matrix of predictors (samples x features), typically from
#'   \code{POSEIDON_encode_for_regression()}.
#' @param method Character string. Scaling method:
#'   \itemize{
#'     \item "zscore" (default): Standardize to mean=0, sd=1
#'     \item "minmax": Scale to 0-1 range
#'     \item "robust": Use median and IQR (robust to outliers)
#'   }
#' @param exclude_binary Logical. If TRUE (default), columns that are binary
#'   (only 0 and 1 values) are not scaled. Dummy-coded categoricals are already
#'   on a meaningful 0/1 scale.
#' @param var_mapping Optional. The var_mapping from \code{POSEIDON_encode_for_regression()}.
#'   If provided, can use variable type information to guide scaling decisions.
#' @param verbose Logical. Print scaling summary (default = TRUE).
#'
#' @return A list containing:
#' \describe{
#'   \item{X_scaled}{Scaled predictor matrix}
#'   \item{scaling_params}{List of parameters for each scaled column:
#'     \itemize{
#'       \item center: Value subtracted (mean, min, or median)
#'       \item scale: Value divided by (sd, range, or IQR)
#'     }
#'     These can be used to transform new data consistently.
#'   }
#'   \item{scaled_cols}{Names of columns that were scaled}
#'   \item{unscaled_cols}{Names of columns left unscaled (binary)}
#'   \item{method}{Scaling method used}
#' }
#'
#' @details
#' Why scaling matters for penalized regression:
#' \itemize{
#'   \item LASSO/Group LASSO apply the same penalty to all coefficients
#'   \item Without scaling, a variable measured in thousands would have a
#'     tiny coefficient (to compensate) and thus be penalized less
#'   \item Scaling puts all variables on equal footing for fair penalization
#' }
#'
#' Binary columns (dummy variables) are typically excluded from scaling because:
#' \itemize{
#'   \item They're already on a 0/1 scale with clear interpretation
#'   \item Scaling would change their interpretation (no longer "presence/absence")
#'   \item For Group LASSO, the group penalty handles their contribution
#' }
#'
#' The scaling parameters are returned so you can apply the same transformation
#' to new/test data, ensuring consistency.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage after encoding
#' encoded <- POSEIDON_encode_for_regression(data, target_var = "outcome")
#' scaled <- POSEIDON_scale_predictors(encoded$X)
#'
#' # Use robust scaling for data with outliers
#' scaled <- POSEIDON_scale_predictors(encoded$X, method = "robust")
#'
#' # Scale everything including binary (not typical)
#' scaled <- POSEIDON_scale_predictors(encoded$X, exclude_binary = FALSE)
#'
#' }
POSEIDON_scale_predictors <- function(X,
                                       method = "zscore",
                                       exclude_binary = TRUE,
                                       var_mapping = NULL,
                                       verbose = TRUE) {

  # Validate inputs
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a matrix or data.frame")
  }

  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }

  method <- tolower(method)
  if (!method %in% c("zscore", "minmax", "robust")) {
    stop("method must be 'zscore', 'minmax', or 'robust'")
  }

  n_samples <- nrow(X)
  n_cols <- ncol(X)
  col_names <- colnames(X)

  if (is.null(col_names)) {
    col_names <- paste0("V", seq_len(n_cols))
    colnames(X) <- col_names
  }

  # Identify binary columns
  is_binary <- apply(X, 2, function(col) {
    unique_vals <- unique(col[!is.na(col)])
    length(unique_vals) <= 2 && all(unique_vals %in% c(0, 1))
  })

  # Determine which columns to scale
  if (exclude_binary) {
    cols_to_scale <- !is_binary
  } else {
    cols_to_scale <- rep(TRUE, n_cols)
  }

  scaled_cols <- col_names[cols_to_scale]
  unscaled_cols <- col_names[!cols_to_scale]

  if (verbose) {
    cat("[POSEIDON] Scaling predictors:\n")
    cat("  Method:", method, "\n")
    cat("  Total columns:", n_cols, "\n")
    cat("  Columns to scale:", sum(cols_to_scale), "\n")
    cat("  Binary columns (unscaled):", sum(!cols_to_scale), "\n")
  }

  # Initialize output
  X_scaled <- X
  scaling_params <- list()

  # Scale each column that needs scaling
  for (j in which(cols_to_scale)) {
    col <- X[, j]
    col_name <- col_names[j]

    # Handle columns with no variance
    col_var <- var(col, na.rm = TRUE)
    if (is.na(col_var) || col_var == 0) {
      if (verbose) {
        warning("Column '", col_name, "' has zero variance, setting to 0")
      }
      X_scaled[, j] <- 0
      scaling_params[[col_name]] <- list(center = mean(col, na.rm = TRUE),
                                          scale = 1,
                                          zero_variance = TRUE)
      next
    }

    if (method == "zscore") {
      center_val <- mean(col, na.rm = TRUE)
      scale_val <- sd(col, na.rm = TRUE)

    } else if (method == "minmax") {
      min_val <- min(col, na.rm = TRUE)
      max_val <- max(col, na.rm = TRUE)
      center_val <- min_val
      scale_val <- max_val - min_val

    } else if (method == "robust") {
      center_val <- median(col, na.rm = TRUE)
      # IQR, but use MAD if IQR is 0
      scale_val <- IQR(col, na.rm = TRUE)
      if (scale_val == 0) {
        scale_val <- mad(col, na.rm = TRUE)
        if (scale_val == 0) scale_val <- 1  # Last resort
      }
    }

    # Apply scaling
    X_scaled[, j] <- (col - center_val) / scale_val

    # Store parameters
    scaling_params[[col_name]] <- list(
      center = center_val,
      scale = scale_val,
      zero_variance = FALSE
    )
  }

  # Add entries for unscaled columns (for completeness)
  for (col_name in unscaled_cols) {
    scaling_params[[col_name]] <- list(
      center = 0,
      scale = 1,
      binary = TRUE
    )
  }

  if (verbose) {
    # Summary statistics of scaled columns
    if (length(scaled_cols) > 0) {
      scaled_means <- colMeans(X_scaled[, scaled_cols, drop = FALSE], na.rm = TRUE)
      scaled_sds <- apply(X_scaled[, scaled_cols, drop = FALSE], 2, sd, na.rm = TRUE)
      cat("\n[POSEIDON] Scaled columns summary:\n")
      cat("  Mean range: [", round(min(scaled_means), 4), ", ",
          round(max(scaled_means), 4), "]\n", sep = "")
      cat("  SD range: [", round(min(scaled_sds), 4), ", ",
          round(max(scaled_sds), 4), "]\n", sep = "")
    }
  }

  return(list(
    X_scaled = X_scaled,
    scaling_params = scaling_params,
    scaled_cols = scaled_cols,
    unscaled_cols = unscaled_cols,
    method = method
  ))
}
