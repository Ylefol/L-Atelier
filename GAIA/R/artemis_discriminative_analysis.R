#' Fit Group LASSO for Variable Selection
#'
#' @description Fits a Group LASSO model for variable selection in mixed data.
#' Group LASSO treats all dummy-coded columns from a single categorical variable
#' as a group, selecting entire variables rather than individual levels.
#'
#' @param X Numeric matrix of predictors (samples x features), typically from
#'   \code{POSEIDON_scale_predictors()}.
#' @param y Response vector. For binary outcomes, should be 0/1 or a factor.
#'   For continuous outcomes, should be numeric.
#' @param groups Integer vector mapping each column in X to its variable group.
#'   Columns with the same group number belong to the same original variable.
#'   Typically from \code{encoded$var_mapping$groups}.
#' @param family Character string. Response type:
#'   \itemize{
#'     \item "binomial": Binary outcome (logistic regression)
#'     \item "multinomial": Multi-class outcome (3+ categories)
#'     \item "gaussian": Continuous outcome (linear regression)
#'     \item "poisson": Count outcome (Poisson regression)
#'   }
#' @param penalty Character string. Penalty type:
#'   \itemize{
#'     \item "grLasso" (default): Group LASSO - L1 penalty on group norms
#'     \item "grMCP": Group MCP - non-convex, less biased estimates
#'     \item "grSCAD": Group SCAD - non-convex alternative
#'   }
#' @param nlambda Integer. Number of lambda values in the regularization path
#'   (default = 100).
#' @param lambda Optional numeric vector. Specific lambda values to use.
#'   If NULL (default), automatically generated.
#' @param alpha Elastic net mixing parameter (default = 1 for pure Group LASSO).
#'   Values < 1 add L2 penalty for stability with correlated predictors.
#' @param group_names Optional character vector. Names for each group (for
#'   nicer output). If NULL, uses "Group1", "Group2", etc.
#' @param verbose Logical. Print fitting summary (default = TRUE).
#'
#' @return An object of class "artemis_group_lasso" containing:
#' \describe{
#'   \item{model}{The fitted grpreg model object}
#'   \item{groups}{Group assignments used}
#'   \item{group_names}{Names for each group}
#'   \item{family}{Response family used}
#'   \item{penalty}{Penalty type used}
#'   \item{n_groups}{Number of variable groups}
#'   \item{n_samples}{Number of samples}
#'   \item{n_features}{Number of features (columns in X)}
#'   \item{lambda}{Vector of lambda values in path}
#' }
#'
#' @details
#' Group LASSO is ideal for variable selection with categorical predictors because:
#' \itemize{
#'   \item It treats all dummy columns from one variable as a unit
#'   \item A variable is either fully in or fully out of the model
#'   \item This avoids the awkward situation where only some levels are selected
#'   \item Results are interpretable: "Variable X is predictive of Y"
#' }
#'
#' The regularization path fits models across a range of lambda values, from
#' high (few/no variables selected) to low (many variables selected). Use
#' \code{ARTEMIS_select_lambda()} to choose the optimal lambda via cross-validation.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # After encoding and scaling
#' encoded <- POSEIDON_encode_for_regression(data, target_var = "outcome")
#' scaled <- POSEIDON_scale_predictors(encoded$X)
#'
#' # Fit Group LASSO for binary outcome
#' fit <- ARTEMIS_fit_group_lasso(
#'   X = scaled$X_scaled,
#'   y = encoded$y,
#'   groups = encoded$var_mapping$groups,
#'   family = "binomial"
#' )
#'
#' # With group names for better output
#' fit <- ARTEMIS_fit_group_lasso(
#'   X = scaled$X_scaled,
#'   y = encoded$y,
#'   groups = encoded$var_mapping$groups,
#'   group_names = encoded$var_mapping$group_names,
#'   family = "binomial"
#' )
#'
#' }
ARTEMIS_fit_group_lasso <- function(X,
                                     y,
                                     groups,
                                     family = "binomial",
                                     penalty = "grLasso",
                                     nlambda = 100,
                                     lambda = NULL,
                                     alpha = 1,
                                     group_names = NULL,
                                     verbose = TRUE) {

  # Check for grpreg package

  if (!requireNamespace("grpreg", quietly = TRUE)) {
    stop("Package 'grpreg' is required for Group LASSO. ",
         "Install with: install.packages('grpreg')")
  }

  # Validate inputs
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

  if (nrow(X) != length(y)) {
    stop("Number of rows in X (", nrow(X), ") must match length of y (",
         length(y), ")")
  }

  if (length(groups) != ncol(X)) {
    stop("Length of groups (", length(groups), ") must match number of columns in X (",
         ncol(X), ")")
  }

  # Detailed validation to catch common grpreg issues
  # Check for NAs in X
  if (any(is.na(X))) {
    na_cols <- which(colSums(is.na(X)) > 0)
    stop("X contains NA values in columns: ",
         paste(head(colnames(X)[na_cols], 5), collapse = ", "),
         if (length(na_cols) > 5) paste0(" (and ", length(na_cols) - 5, " more)") else "",
         "\ngrpreg cannot handle missing values. Remove or impute NAs first.")
  }

  # Check for NAs in y
  if (any(is.na(y))) {
    stop("y contains ", sum(is.na(y)), " NA values. ",
         "Remove rows with NA in target before fitting.")
  }

  # Check for infinite values
  if (any(is.infinite(X))) {
    inf_cols <- which(colSums(is.infinite(X)) > 0)
    stop("X contains infinite values in columns: ",
         paste(head(colnames(X)[inf_cols], 5), collapse = ", "))
  }

  # Check that X is numeric
  if (!is.numeric(X)) {
    stop("X must be numeric. Current storage mode: ", storage.mode(X))
  }

  # Check that groups are valid integers
  if (!is.numeric(groups)) {
    stop("groups must be numeric. Current class: ", class(groups)[1])
  }

  # Ensure groups are consecutive integers starting from 1
  unique_groups <- sort(unique(groups))
  expected_groups <- seq_len(max(groups))
  if (!all(unique_groups == expected_groups)) {
    missing_groups <- setdiff(expected_groups, unique_groups)
    stop("groups must be consecutive integers from 1 to max. ",
         "Missing group numbers: ", paste(head(missing_groups, 5), collapse = ", "))
  }

  # Check for constant columns (zero variance)
  col_vars <- apply(X, 2, var, na.rm = TRUE)
  zero_var_cols <- which(col_vars == 0 | is.na(col_vars))
  if (length(zero_var_cols) > 0) {
    warning("X contains ", length(zero_var_cols), " constant column(s) (zero variance): ",
            paste(head(colnames(X)[zero_var_cols], 5), collapse = ", "),
            if (length(zero_var_cols) > 5) paste0(" (and ", length(zero_var_cols) - 5, " more)") else "",
            "\nThese may cause issues. Consider removing them.")
  }

  # Check y is not constant

  if (length(unique(y)) < 2) {
    stop("y has only ", length(unique(y)), " unique value(s). ",
         "Need at least 2 distinct values for regression.")
  }

  # Validate family
  family <- tolower(family)
  if (!family %in% c("binomial", "multinomial", "gaussian", "poisson")) {
    stop("family must be 'binomial', 'multinomial', 'gaussian', or 'poisson'")
  }

  # Validate penalty
  penalty <- tolower(penalty)
  # grpreg uses specific capitalization
  penalty_map <- c("grlasso" = "grLasso", "grmcp" = "grMCP", "grscad" = "grSCAD")
  if (!penalty %in% names(penalty_map)) {
    stop("penalty must be 'grLasso', 'grMCP', or 'grSCAD'")
  }
  penalty <- penalty_map[penalty]

  # Handle response variable
  if (family == "binomial") {
    if (is.factor(y)) {
      y <- as.numeric(y) - 1  # Convert factor to 0/1
    }
    if (!all(y %in% c(0, 1))) {
      stop("For family='binomial', y must be binary (0/1 or a 2-level factor)")
    }
  } else if (family == "multinomial") {
    if (!is.factor(y)) {
      y <- as.factor(y)
    }
    n_classes <- nlevels(y)
    if (n_classes < 3) {
      warning("Only ", n_classes, " classes detected. Consider using 'binomial' for binary outcomes.")
    }
    if (verbose) {
      cat("[ARTEMIS] Classes:", paste(levels(y), collapse = ", "), "\n")
    }
  }

  # Set up group names
  n_groups <- max(groups)
  if (is.null(group_names)) {
    group_names <- paste0("Group", seq_len(n_groups))
  } else if (length(group_names) != n_groups) {
    warning("Length of group_names doesn't match number of groups. Using defaults.")
    group_names <- paste0("Group", seq_len(n_groups))
  }

  n_samples <- nrow(X)
  n_features <- ncol(X)

  if (verbose) {
    cat("[ARTEMIS] Fitting Group LASSO:\n")
    cat("    Samples:", n_samples, "\n")
    cat("    Features:", n_features, "\n")
    cat("    Groups:", n_groups, "\n")
    cat("    Family:", family, "\n")
    cat("    Penalty:", penalty, "\n")
    cat("    Lambda values:", nlambda, "\n")
  }

  # Fit the model
  if (verbose) cat("[ARTEMIS] Fitting regularization path...\n")

  # grpreg can be sensitive to named vectors and attributes
  # Strip names and ensure clean vectors
  X_clean <- unname(X)
  # For multinomial, keep as factor; for others, convert to numeric
  if (family == "multinomial") {
    y_clean <- factor(y)
  } else {
    y_clean <- as.numeric(unname(y))
  }
  groups_clean <- as.integer(unname(groups))

  # Debug output if verbose
  if (verbose) {
    cat("    Data check:\n")
    cat("        X: ", nrow(X_clean), "x", ncol(X_clean), ", class=", class(X_clean)[1],
        ", mode=", storage.mode(X_clean), "\n", sep = "")
    if (family == "multinomial") {
      cat("        y: length=", length(y_clean), ", classes=", nlevels(y_clean), "\n", sep = "")
    } else {
      cat("        y: length=", length(y_clean), ", class=", class(y_clean)[1],
          ", unique values=", length(unique(y_clean)), "\n", sep = "")
    }
    cat("        groups: length=", length(groups_clean), ", range=[", min(groups_clean),
        ",", max(groups_clean), "]\n", sep = "")
  }

  # Build arguments list - only include non-NULL optional parameters
  # grpreg can fail with explicit NULL arguments
  grpreg_args <- list(
    X = X_clean,
    y = y_clean,
    group = groups_clean,
    penalty = penalty,
    family = family,
    nlambda = nlambda
  )

  # Only add lambda if user specified it
  if (!is.null(lambda)) {
    grpreg_args$lambda <- lambda
  }

  # Only add alpha if not default (1)
  if (alpha != 1) {
    grpreg_args$alpha <- alpha
  }

  fit <- do.call(grpreg::grpreg, grpreg_args)

  if (verbose) {
    cat("    Done.\n")
    cat("    Lambda range: [", round(min(fit$lambda), 4), ", ",
        round(max(fit$lambda), 4), "]\n", sep = "")

    # Count variables selected at different points
    # Use coef() directly - it's a base R generic that grpreg extends
    n_selected_min <- sum(coef(fit, which = length(fit$lambda))[-1] != 0)
    n_selected_max <- sum(coef(fit, which = 1)[-1] != 0)
    cat("    Features selected: ", n_selected_max, " (lambda.max) to ",
        n_selected_min, " (lambda.min)\n", sep = "")
  }

  # Build result object
  result <- list(
    model = fit,
    groups = groups,
    group_names = group_names,
    family = family,
    penalty = penalty,
    n_groups = n_groups,
    n_samples = n_samples,
    n_features = n_features,
    lambda = fit$lambda
  )

  class(result) <- c("artemis_group_lasso", "list")

  return(result)
}


#' Print Method for artemis_group_lasso
#'
#' @param x An artemis_group_lasso object
#' @param ... Additional arguments (unused)
#' @method print artemis_group_lasso
#' @export
print.artemis_group_lasso <- function(x, ...) {
  cat("Group LASSO Model (artemis_group_lasso)\n")
  cat("------------------------------\n")
  cat("Samples:", x$n_samples, "\n")
  cat("Features:", x$n_features, "\n")
  cat("Groups:", x$n_groups, "\n")
  cat("Family:", x$family, "\n")
  cat("Penalty:", x$penalty, "\n")
  cat("Lambda values:", length(x$lambda), "\n")
  cat("Lambda range: [", round(min(x$lambda), 4), ", ",
      round(max(x$lambda), 4), "]\n", sep = "")
  cat("\nUse ARTEMIS_select_lambda() to choose optimal lambda via CV\n")
  invisible(x)
}


#' Select Optimal Lambda via Cross-Validation
#'
#' @description Performs cross-validation to select the optimal regularization
#' parameter (lambda) for a Group LASSO model. Returns both lambda.min (minimum
#' CV error) and lambda.1se (most parsimonious within 1 SE of minimum).
#'
#' @param glasso_fit An artemis_group_lasso object from \code{ARTEMIS_fit_group_lasso()}.
#' @param X The predictor matrix used to fit the model.
#' @param y The response vector used to fit the model.
#' @param nfolds Integer. Number of cross-validation folds (default = 10).
#' @param seed Integer. Random seed for reproducibility (default = NULL).
#' @param verbose Logical. Print CV summary (default = TRUE).
#'
#' @return A list containing:
#' \describe{
#'   \item{lambda.min}{Lambda with minimum CV error}
#'   \item{lambda.1se}{Largest lambda within 1 SE of minimum (more parsimonious)}
#'   \item{cve}{Cross-validation error for each lambda}
#'   \item{cvse}{Standard error of CV error}
#'   \item{cv_fit}{The cv.grpreg object for further inspection}
#'   \item{n_selected_min}{Number of groups selected at lambda.min}
#'   \item{n_selected_1se}{Number of groups selected at lambda.1se}
#' }
#'
#' @details
#' Two lambda choices are provided:
#' \itemize{
#'   \item \code{lambda.min}: Minimizes prediction error. May include more
#'     variables than necessary.
#'   \item \code{lambda.1se}: The largest lambda (most regularized) whose error
#'     is within 1 standard error of the minimum. More parsimonious, recommended
#'     for variable selection.
#' }
#'
#' For variable selection purposes, \code{lambda.1se} is typically preferred as
#' it provides a simpler model with similar predictive performance.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # After fitting Group LASSO
#' fit <- ARTEMIS_fit_group_lasso(X, y, groups, family = "binomial")
#'
#' # Select optimal lambda
#' cv_result <- ARTEMIS_select_lambda(fit, X, y, nfolds = 10)
#'
#' # Use lambda.1se for variable selection (recommended)
#' selected <- ARTEMIS_extract_selected_variables(fit, cv_result$lambda.1se)
#'
#' }
ARTEMIS_select_lambda <- function(glasso_fit,
                                   X,
                                   y,
                                   nfolds = 10,
                                   seed = NULL,
                                   verbose = TRUE) {

  # Validate input

  if (!inherits(glasso_fit, "artemis_group_lasso")) {
    stop("glasso_fit must be an artemis_group_lasso object from ARTEMIS_fit_group_lasso()")
  }

  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

  # Handle response based on family
  if (glasso_fit$family == "binomial" && is.factor(y)) {
    y <- as.numeric(y) - 1
  } else if (glasso_fit$family == "multinomial" && !is.factor(y)) {
    y <- as.factor(y)
  }

  if (verbose) {
    cat("[ARTEMIS] Cross-validation for lambda selection:\n")
    cat("    Folds:", nfolds, "\n")
  }

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Run cross-validation
  cv_fit <- grpreg::cv.grpreg(
    X = X,
    y = y,
    group = glasso_fit$groups,
    penalty = glasso_fit$penalty,
    family = glasso_fit$family,
    nfolds = nfolds,
    lambda = glasso_fit$lambda
  )

  # Extract results
  lambda_min <- cv_fit$lambda.min

  # Calculate lambda.1se manually if not provided
  # (some versions of grpreg don't have it)
  cve <- cv_fit$cve
  cvse <- cv_fit$cvse
  min_idx <- which.min(cve)
  threshold <- cve[min_idx] + cvse[min_idx]

  # Find largest lambda (smallest index) with error <= threshold
  candidates <- which(cve <= threshold)
  lambda_1se_idx <- min(candidates)
  lambda_1se <- cv_fit$lambda[lambda_1se_idx]

  # Count selected groups at each lambda
  # Use coef() directly - it's a base R generic that grpreg extends
  coef_min <- coef(glasso_fit$model, lambda = lambda_min)[-1]  # Remove intercept
  coef_1se <- coef(glasso_fit$model, lambda = lambda_1se)[-1]

  # Count groups with at least one non-zero coefficient
  groups <- glasso_fit$groups
  n_groups <- glasso_fit$n_groups

  count_selected_groups <- function(coefs, groups) {
    selected <- 0
    for (g in seq_len(max(groups))) {
      group_coefs <- coefs[groups == g]
      if (any(group_coefs != 0)) {
        selected <- selected + 1
      }
    }
    return(selected)
  }

  n_selected_min <- count_selected_groups(coef_min, groups)
  n_selected_1se <- count_selected_groups(coef_1se, groups)

  if (verbose) {
    cat("    Lambda.min:", round(lambda_min, 4),
        "(", n_selected_min, "groups selected)\n")
    cat("    Lambda.1se:", round(lambda_1se, 4),
        "(", n_selected_1se, "groups selected)\n")
    cat("\n    Recommendation: Use lambda.1se for more parsimonious selection\n")
  }

  return(list(
    lambda.min = lambda_min,
    lambda.1se = lambda_1se,
    cve = cve,
    cvse = cvse,
    cv_fit = cv_fit,
    n_selected_min = n_selected_min,
    n_selected_1se = n_selected_1se
  ))
}


#' Extract Selected Variables from Group LASSO
#'
#' @description Identifies which variable groups have non-zero coefficients
#' at a given lambda value. Returns the names of selected variables.
#'
#' @param glasso_fit An artemis_group_lasso object from \code{ARTEMIS_fit_group_lasso()}.
#' @param lambda Numeric. The lambda value at which to extract selected variables.
#'   Typically \code{cv_result$lambda.1se} or \code{cv_result$lambda.min}.
#' @param verbose Logical. Print selection summary (default = TRUE).
#'
#' @return A list containing:
#' \describe{
#'   \item{selected_groups}{Integer vector of selected group numbers}
#'   \item{selected_names}{Character vector of selected variable names}
#'   \item{n_selected}{Number of variables selected}
#'   \item{n_total}{Total number of variables}
#'   \item{lambda}{Lambda value used}
#'   \item{coefficients}{Named vector of non-zero coefficients (all columns)}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # After fitting and CV
#' fit <- ARTEMIS_fit_group_lasso(X, y, groups, family = "binomial")
#' cv_result <- ARTEMIS_select_lambda(fit, X, y)
#'
#' # Extract selected variables at lambda.1se
#' selected <- ARTEMIS_extract_selected_variables(fit, cv_result$lambda.1se)
#' print(selected$selected_names)
#'
#' }
ARTEMIS_extract_selected_variables <- function(glasso_fit,
                                                lambda,
                                                verbose = TRUE) {

  if (!inherits(glasso_fit, "artemis_group_lasso")) {
    stop("glasso_fit must be an artemis_group_lasso object")
  }

  # Get coefficients at specified lambda
  coefs_raw <- coef(glasso_fit$model, lambda = lambda)

  # Handle multinomial vs other families
  # Multinomial returns a matrix (features x classes), others return a vector
  if (glasso_fit$family == "multinomial") {
    # For multinomial, coefs_raw is a matrix with intercepts in first row
    intercept <- coefs_raw[1, ]
    coefs_matrix <- coefs_raw[-1, , drop = FALSE]
    # A variable is selected if ANY class has non-zero coefficient
    coefs <- apply(abs(coefs_matrix), 1, max)  # Max absolute coef across classes
    nonzero_mask <- coefs != 0
  } else {
    # For binomial/gaussian/poisson, coefs_raw is a vector
    intercept <- coefs_raw[1]
    coefs <- coefs_raw[-1]
    nonzero_mask <- coefs != 0
  }

  # Identify which groups have at least one non-zero coefficient
  groups <- glasso_fit$groups
  group_names <- glasso_fit$group_names
  n_groups <- glasso_fit$n_groups

  selected_groups <- integer(0)
  for (g in seq_len(n_groups)) {
    group_mask <- groups == g
    if (any(nonzero_mask[group_mask])) {
      selected_groups <- c(selected_groups, g)
    }
  }

  selected_names <- group_names[selected_groups]
  n_selected <- length(selected_groups)

  # Get non-zero coefficients for reporting
  nonzero_coefs <- coefs[nonzero_mask]
  names(nonzero_coefs) <- names(coefs)[nonzero_mask]

  if (verbose) {
    cat("[ARTEMIS] Selected variables at lambda =", round(lambda, 4), ":\n")
    cat("    ", n_selected, " of ", n_groups, " variables selected\n\n", sep = "")

    if (n_selected > 0) {
      cat("    Selected variables:\n")
      for (i in seq_along(selected_names)) {
        var_name <- selected_names[i]
        group_idx <- selected_groups[i]
        group_coefs <- coefs[groups == group_idx]
        n_nonzero <- sum(group_coefs != 0)
        n_total_cols <- length(group_coefs)

        if (n_total_cols == 1) {
          # Quantitative variable - show coefficient
          cat("    ", i, ". ", var_name, " (coef = ",
              round(group_coefs[group_coefs != 0], 4), ")\n", sep = "")
        } else {
          # Categorical variable - show how many levels have non-zero coefs
          cat("    ", i, ". ", var_name, " (", n_nonzero, "/", n_total_cols,
              " levels active)\n", sep = "")
        }
      }
    }
  }

  return(list(
    selected_groups = selected_groups,
    selected_names = selected_names,
    n_selected = n_selected,
    n_total = n_groups,
    lambda = lambda,
    intercept = intercept,
    coefficients = nonzero_coefs
  ))
}


#' Fit Final Interpretable Model on Selected Variables
#'
#' @description After variable selection with Group LASSO, fits an unpenalized
#' regression model on only the selected variables. This provides unbiased
#' coefficient estimates with proper p-values and confidence intervals for
#' inference and interpretation.
#'
#' @param data Original data.frame (before encoding). Must contain the target
#'   variable and all selected predictor variables.
#' @param target_var Character string. Name of the target variable column.
#' @param selected_vars Character vector. Names of selected variables from
#'   \code{ARTEMIS_extract_selected_variables()$selected_names}.
#' @param family Character string. Response type:
#'   \itemize{
#'     \item "binomial": Binary outcome (logistic regression)
#'     \item "multinomial": Multi-class outcome (3+ categories)
#'     \item "gaussian": Continuous outcome (linear regression)
#'   }
#' @param verbose Logical. Print model summary (default = TRUE).
#'
#' @return A list of class "artemis_final_model" containing:
#' \describe{
#'   \item{model}{The fitted glm object}
#'   \item{formula}{The model formula used}
#'   \item{family}{Response family}
#'   \item{selected_vars}{Variables included in model}
#'   \item{n_samples}{Number of samples used}
#'   \item{converged}{Whether the model converged}
#' }
#'
#' @details
#' Why refit without penalty?
#' \itemize{
#'   \item Group LASSO coefficients are shrunk toward zero (biased)
#'   \item P-values and CIs from penalized models are not valid for inference
#'   \item Refitting on selected variables gives unbiased estimates
#'   \item Standard errors, p-values, and CIs are now interpretable
#' }
#'
#' This is the "post-selection inference" step. The selected variables were
#' chosen by Group LASSO; now we estimate their effects properly.
#'
#' Note: P-values should be interpreted cautiously since variables were
#' pre-selected. The selection step already established these variables are
#' predictive; the final model quantifies their effects.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # After selection
#' selected <- ARTEMIS_extract_selected_variables(fit, cv_result$lambda.1se)
#'
#' # Fit final model on original data
#' final <- ARTEMIS_fit_final_model(
#'   data = clinical_data,
#'   target_var = "outcome",
#'   selected_vars = selected$selected_names,
#'   family = "binomial"
#' )
#'
#' # View standard model summary
#' summary(final$model)
#'
#' }
ARTEMIS_fit_final_model <- function(data,
                                     target_var,
                                     selected_vars,
                                     family = "binomial",
                                     verbose = TRUE) {

  # Validate inputs
  if (!is.data.frame(data)) {
    stop("data must be a data.frame")
  }

  if (!target_var %in% colnames(data)) {
    stop("target_var '", target_var, "' not found in data")
  }

  missing_vars <- setdiff(selected_vars, colnames(data))
  if (length(missing_vars) > 0) {
    stop("Selected variables not found in data: ",
         paste(missing_vars, collapse = ", "))
  }

  if (length(selected_vars) == 0) {
    stop("No variables selected. Cannot fit model with no predictors.")
  }

  family <- tolower(family)
  if (!family %in% c("binomial", "multinomial", "gaussian")) {
    stop("family must be 'binomial', 'multinomial', or 'gaussian'")
  }

  # Check for nnet package if multinomial
  if (family == "multinomial") {
    if (!requireNamespace("nnet", quietly = TRUE)) {
      stop("Package 'nnet' is required for multinomial regression. ",
           "Install with: install.packages('nnet')")
    }
  }

  # Build formula
  predictor_string <- paste(selected_vars, collapse = " + ")
  formula_string <- paste(target_var, "~", predictor_string)
  model_formula <- as.formula(formula_string)

  if (verbose) {
    cat("[ARTEMIS] Fitting final model:\n")
    cat("    Target:", target_var, "\n")
    cat("    Predictors:", length(selected_vars), "variables\n")
    cat("    Family:", family, "\n")
    cat("    Formula:", formula_string, "\n\n")
  }

  # Subset data to relevant columns (avoid issues with other columns)
  model_data <- data[, c(target_var, selected_vars), drop = FALSE]

  # Remove rows with NA in any of the model variables
  complete_rows <- complete.cases(model_data)
  n_removed <- sum(!complete_rows)
  if (n_removed > 0) {
    if (verbose) {
      cat("    Removed", n_removed, "rows with missing values\n")
    }
    model_data <- model_data[complete_rows, , drop = FALSE]
  }

  n_samples <- nrow(model_data)

  if (verbose) {
    cat("    Samples:", n_samples, "\n\n")
  }

  # Fit the model
  if (family == "binomial") {
    model <- glm(model_formula, data = model_data, family = binomial(link = "logit"))
    converged <- model$converged
  } else if (family == "multinomial") {
    # Ensure target is a factor
    model_data[[target_var]] <- as.factor(model_data[[target_var]])
    # Use nnet::multinom for multinomial logistic regression
    # trace = FALSE suppresses iteration output
    model <- nnet::multinom(model_formula, data = model_data, trace = FALSE)
    converged <- model$convergence == 0
    if (verbose) {
      cat("    Reference class:", levels(model_data[[target_var]])[1], "\n")
    }
  } else {
    model <- glm(model_formula, data = model_data, family = gaussian(link = "identity"))
    converged <- model$converged
  }

  if (!converged && verbose) {
    warning("Model did not converge. Results may be unreliable.")
  }

  if (verbose) {
    cat("[ARTEMIS] Model fit complete.\n")
    if (family == "binomial") {
      cat("    Null deviance:", round(model$null.deviance, 2), "on",
          model$df.null, "df\n")
      cat("    Residual deviance:", round(model$deviance, 2), "on",
          model$df.residual, "df\n")
      cat("    AIC:", round(model$aic, 2), "\n")
    } else if (family == "multinomial") {
      cat("    Residual deviance:", round(model$deviance, 2), "\n")
      cat("    AIC:", round(model$AIC, 2), "\n")
      cat("    Classes:", paste(model$lev, collapse = ", "), "\n")
    } else if (family == "gaussian") {
      # For glm with gaussian family, sigma is sqrt of dispersion
      model_summary <- summary(model)
      residual_se <- sqrt(model_summary$dispersion)
      cat("    Residual SE:", round(residual_se, 4), "\n")
      cat("    R-squared:", round(1 - model$deviance/model$null.deviance, 4), "\n")
    } else {
      # Poisson or other families
      cat("    Residual deviance:", round(model$deviance, 2), "on",
          model$df.residual, "df\n")
      cat("    AIC:", round(model$aic, 2), "\n")
    }
  }

  result <- list(
    model = model,
    formula = model_formula,
    family = family,
    selected_vars = selected_vars,
    n_samples = n_samples,
    converged = converged
  )

  class(result) <- c("artemis_final_model", "list")

  return(result)
}


#' Extract Coefficients with Confidence Intervals
#'
#' @description Extracts model coefficients in an interpretable format with
#' confidence intervals and p-values. For logistic regression, converts
#' coefficients to odds ratios. Returns a publication-ready table.
#'
#' @param final_model An artemis_final_model object from \code{ARTEMIS_fit_final_model()},
#'   or a standard glm object.
#' @param format Character string. Output format:
#'   \itemize{
#'     \item "odds_ratio" (default for binomial): Exponentiated coefficients
#'     \item "coefficient": Raw coefficients (log-odds for binomial)
#'   }
#' @param conf_level Numeric. Confidence level for intervals (default = 0.95).
#' @param digits Integer. Number of decimal places for rounding (default = 3).
#' @param include_intercept Logical. Include intercept in output (default = FALSE).
#' @param verbose Logical. Print formatted table (default = TRUE).
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{variable}{Variable name (and level for categorical)}
#'   \item{estimate}{Coefficient or odds ratio}
#'   \item{ci_lower}{Lower confidence bound}
#'   \item{ci_upper}{Upper confidence bound}
#'   \item{std_error}{Standard error}
#'   \item{z_value}{Z-statistic (or t for gaussian)}
#'   \item{p_value}{P-value}
#'   \item{significance}{Significance stars (* p<0.05, ** p<0.01, *** p<0.001)}
#' }
#'
#' @details
#' For logistic regression (family = binomial):
#' \itemize{
#'   \item Odds ratio > 1: Higher values of predictor associated with higher
#'     probability of outcome = 1
#'   \item Odds ratio < 1: Higher values associated with lower probability
#'   \item Odds ratio = 1: No association
#' }
#'
#' For categorical variables, each level (except reference) gets its own row.
#' The odds ratio compares that level to the reference level.
#'
#' Confidence intervals are Wald-based (coefficient ± z * SE).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # After fitting final model
#' final <- ARTEMIS_fit_final_model(data, "outcome", selected_vars, "binomial")
#'
#' # Get odds ratios with 95% CI
#' coef_table <- ARTEMIS_extract_coefficients(final)
#'
#' # Get raw coefficients instead
#' coef_table <- ARTEMIS_extract_coefficients(final, format = "coefficient")
#'
#' # Use 99% confidence intervals
#' coef_table <- ARTEMIS_extract_coefficients(final, conf_level = 0.99)
#'
#' }
ARTEMIS_extract_coefficients <- function(final_model,
                                          format = NULL,
                                          conf_level = 0.95,
                                          digits = 3,
                                          include_intercept = FALSE,
                                          verbose = TRUE) {

  # Handle both artemis_final_model and raw glm/multinom objects
  if (inherits(final_model, "artemis_final_model")) {
    model <- final_model$model
    family <- final_model$family
  } else if (inherits(final_model, "glm")) {
    model <- final_model
    family <- model$family$family
  } else if (inherits(final_model, "multinom")) {
    model <- final_model
    family <- "multinomial"
  } else {
    stop("final_model must be an artemis_final_model, glm, or multinom object")
  }

  # Set default format based on family
  if (is.null(format)) {
    format <- if (family %in% c("binomial", "multinomial")) "odds_ratio" else "coefficient"
  }

  format <- tolower(format)
  if (!format %in% c("odds_ratio", "coefficient")) {
    stop("format must be 'odds_ratio' or 'coefficient'")
  }

  # Handle multinomial separately due to different coefficient structure
  if (family == "multinomial") {
    return(.extract_multinomial_coefficients(model, format, conf_level, digits,
                                              include_intercept, verbose))
  }

  # Get model summary (for binomial/gaussian)
  model_summary <- summary(model)
  coef_table <- as.data.frame(model_summary$coefficients)

  # Standardize column names (different between gaussian and binomial)
  if (ncol(coef_table) == 4) {
    colnames(coef_table) <- c("estimate", "std_error", "statistic", "p_value")
  }

  # Add variable names
  coef_table$variable <- rownames(coef_table)
  rownames(coef_table) <- NULL

  # Remove intercept if requested
  if (!include_intercept) {
    coef_table <- coef_table[coef_table$variable != "(Intercept)", , drop = FALSE]
  }

  # Calculate confidence intervals
  z_value <- qnorm(1 - (1 - conf_level) / 2)

  coef_table$ci_lower <- coef_table$estimate - z_value * coef_table$std_error
  coef_table$ci_upper <- coef_table$estimate + z_value * coef_table$std_error

  # Transform to odds ratios if requested
  if (format == "odds_ratio") {
    coef_table$estimate <- exp(coef_table$estimate)
    coef_table$ci_lower <- exp(coef_table$ci_lower)
    coef_table$ci_upper <- exp(coef_table$ci_upper)
  }

  # Add significance stars
  coef_table$significance <- ""
  coef_table$significance[coef_table$p_value < 0.05] <- "*"
  coef_table$significance[coef_table$p_value < 0.01] <- "**"
  coef_table$significance[coef_table$p_value < 0.001] <- "***"

  # Round values
  numeric_cols <- c("estimate", "ci_lower", "ci_upper", "std_error", "statistic")
  for (col in numeric_cols) {
    coef_table[[col]] <- round(coef_table[[col]], digits)
  }

  # Format p-values (keep more precision for small values)
  coef_table$p_value <- ifelse(
    coef_table$p_value < 0.001,
    formatC(coef_table$p_value, format = "e", digits = 2),
    round(coef_table$p_value, digits)
  )

  # Reorder columns
  col_order <- c("variable", "estimate", "ci_lower", "ci_upper",
                 "std_error", "statistic", "p_value", "significance")
  coef_table <- coef_table[, col_order]

  # Rename statistic column based on family
  stat_name <- if (family == "binomial") "z_value" else "t_value"
  colnames(coef_table)[colnames(coef_table) == "statistic"] <- stat_name

  if (verbose) {
    # Print formatted output
    estimate_label <- if (format == "odds_ratio") "Odds Ratio" else "Coefficient"
    ci_label <- paste0(conf_level * 100, "% CI")

    cat("\n")
    cat("[ARTEMIS] Model Coefficients")
    if (format == "odds_ratio") cat("    (Odds Ratios)")
    cat("\n")
    cat("    ", sprintf("%-30s %10s %18s %10s %5s\n",
                "Variable", estimate_label, ci_label, "P-value", ""))
    cat("    ", paste(rep("-", 70), collapse = ""), "\n")

    for (i in seq_len(nrow(coef_table))) {
      row <- coef_table[i, ]
      ci_str <- sprintf("[%s, %s]", row$ci_lower, row$ci_upper)
      cat("    ", sprintf("%-30s %10s %18s %10s %s\n",
                  substr(row$variable, 1, 30),
                  row$estimate,
                  ci_str,
                  row$p_value,
                  row$significance))
    }

    cat("    ", paste(rep("-", 70), collapse = ""), "\n")
    cat("    Signif. codes: *** p<0.001, ** p<0.01, * p<0.05\n\n")
  }

  # Return invisibly if verbose, visibly if not
  if (verbose) {
    invisible(coef_table)
  } else {
    return(coef_table)
  }
}
#' Print method for artemis_final_model
#'
#' @param x An artemis_final_model object
#' @param ... Additional arguments (unused)
#' @method print artemis_final_model
#' @export
print.artemis_final_model <- function(x, ...) {
  cat("Final Model (artemis_final_model)\n")
  cat("------------------------------\n")
  cat("Family:", x$family, "\n")
  cat("Samples:", x$n_samples, "\n")
  cat("Variables:", length(x$selected_vars), "\n")
  cat("Converged:", x$converged, "\n")
  cat("\nFormula:", deparse(x$formula), "\n")
  cat("\nUse summary(obj$model) for coefficient details\n")
  cat("Use ARTEMIS_extract_coefficients() for formatted output\n")
  invisible(x)
}


#' Extract Coefficients from Multinomial Model (Internal)
#'
#' @description Internal helper function to extract and format coefficients
#' from a nnet::multinom model. Returns a table with one row per
#' variable-class combination.
#'
#' @param model A multinom model object
#' @param format "odds_ratio" or "coefficient"
#' @param conf_level Confidence level
#' @param digits Rounding digits
#' @param include_intercept Include intercepts?
#' @param verbose Print output?
#'
#' @return data.frame of coefficients
#' @keywords internal
.extract_multinomial_coefficients <- function(model, format, conf_level, digits,
                                               include_intercept, verbose) {

  # Get coefficients and standard errors
  coefs <- coef(model)
  se <- summary(model)$standard.errors

  # Handle single class vs multiple classes
  # (coef returns matrix for >2 classes, may return vector for 2)
  if (is.null(dim(coefs))) {
    # Convert to matrix format
    coefs <- matrix(coefs, nrow = 1)
    se <- matrix(se, nrow = 1)
    rownames(coefs) <- rownames(se) <- model$lev[2]  # Non-reference class
    colnames(coefs) <- colnames(se) <- names(coef(model))
  }

  class_names <- rownames(coefs)
  var_names <- colnames(coefs)

  # Build results table
  results <- data.frame(
    class = character(),
    variable = character(),
    estimate = numeric(),
    std_error = numeric(),
    ci_lower = numeric(),
    ci_upper = numeric(),
    z_value = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )

  z_crit <- qnorm(1 - (1 - conf_level) / 2)

  for (cls in class_names) {
    for (var in var_names) {
      # Skip intercept if requested
      if (!include_intercept && var == "(Intercept)") next

      est <- coefs[cls, var]
      std_err <- se[cls, var]
      z_val <- est / std_err
      p_val <- 2 * pnorm(-abs(z_val))

      ci_low <- est - z_crit * std_err
      ci_high <- est + z_crit * std_err

      results <- rbind(results, data.frame(
        class = cls,
        variable = var,
        estimate = est,
        std_error = std_err,
        ci_lower = ci_low,
        ci_upper = ci_high,
        z_value = z_val,
        p_value = p_val,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Transform to odds ratios if requested
  if (format == "odds_ratio") {
    results$estimate <- exp(results$estimate)
    results$ci_lower <- exp(results$ci_lower)
    results$ci_upper <- exp(results$ci_upper)
  }

  # Add significance stars
  results$significance <- ""
  results$significance[results$p_value < 0.05] <- "*"
  results$significance[results$p_value < 0.01] <- "**"
  results$significance[results$p_value < 0.001] <- "***"

  # Round values
  results$estimate <- round(results$estimate, digits)
  results$ci_lower <- round(results$ci_lower, digits)
  results$ci_upper <- round(results$ci_upper, digits)
  results$std_error <- round(results$std_error, digits)
  results$z_value <- round(results$z_value, digits)

  # Format p-values
  results$p_value <- ifelse(
    results$p_value < 0.001,
    formatC(results$p_value, format = "e", digits = 2),
    round(results$p_value, digits)
  )

  if (verbose) {
    estimate_label <- if (format == "odds_ratio") "Odds Ratio" else "Coefficient"
    ci_label <- paste0(conf_level * 100, "% CI")

    cat("\n")
    cat("[ARTEMIS] Multinomial Model Coefficients")
    if (format == "odds_ratio") cat("    (Odds Ratios)")
    cat("\n")
    cat("    Reference class:", model$lev[1], "\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")

    for (cls in class_names) {
      cat("    Class:", cls, "vs", model$lev[1], "(reference)\n")
      cat(paste(rep("-", 80), collapse = ""), "\n")
      cat("    ", sprintf("%-25s %10s %18s %10s %5s\n",
                  "Variable", estimate_label, ci_label, "P-value", ""))
      cat(paste(rep("-", 80), collapse = ""), "\n")

      cls_results <- results[results$class == cls, ]
      for (i in seq_len(nrow(cls_results))) {
        row <- cls_results[i, ]
        ci_str <- sprintf("[%s, %s]", row$ci_lower, row$ci_upper)
        cat("    ", sprintf("%-25s %10s %18s %10s %s\n",
                    substr(row$variable, 1, 25),
                    row$estimate,
                    ci_str,
                    row$p_value,
                    row$significance))
      }
      cat("\n")
    }

    cat(paste(rep("-", 80), collapse = ""), "\n")
    cat("    Signif. codes: *** p<0.001, ** p<0.01, * p<0.05\n\n")
  }

  if (verbose) {
    invisible(results)
  } else {
    return(results)
  }
}
