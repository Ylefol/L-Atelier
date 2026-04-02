###############################################################################
########### Discriminative Variable Analysis Report ###########
###############################################################################

#' Generate Discriminative Analysis Report
#'
#' @description Generate a comprehensive report from the discriminative variable
#' analysis pipeline, including variable selection results, model coefficients,
#' and optionally save tables and plots.
#'
#' @param glasso_fit An artemis_group_lasso object from ARTEMIS_fit_group_lasso().
#' @param cv_result Output from ARTEMIS_select_lambda().
#' @param selected Output from ARTEMIS_extract_selected_variables().
#' @param final_model Optional. An artemis_final_model from ARTEMIS_fit_final_model().
#' @param coef_table Optional. Data.frame from ARTEMIS_extract_coefficients().
#' @param encoded Optional. Output from POSEIDON_encode_for_regression() for
#'   additional context.
#' @param output_dir Character. Directory to save output files. If NULL (default),
#'   only prints to console without saving files.
#' @param prefix Character. Prefix for output filenames. Default = "discriminative".
#' @param save_tables Logical. Save tables as CSV files. Default = TRUE.
#' @param save_plots Logical. Save plots as PNG files. Default = TRUE.
#' @param plot_width Numeric. Width of saved plots in inches. Default = 10.
#' @param plot_height Numeric. Height of saved plots in inches. Default = 8.
#' @param verbose Logical. Print report to console. Default = TRUE.
#'
#' @return Invisibly returns a list containing:
#' \describe{
#'   \item{summary}{Character vector of the text report}
#'   \item{tables}{List of data frames exported}
#'   \item{files_saved}{Character vector of files saved}
#' }
#'
#' @details
#' This function consolidates results from the discriminative analysis pipeline:
#' \itemize{
#'   \item Group LASSO variable selection results
#'   \item Cross-validation lambda selection
#'   \item Selected variables and their effects
#'   \item Final model coefficients (odds ratios or coefficients)
#' }
#'
#' At minimum, glasso_fit, cv_result, and selected are required. Final model
#' and coefficients are optional but recommended for complete reporting.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Console report only
#' ELEUTHIA_discriminative_report(glasso_fit, cv_result, selected)
#'
#' # Full report with file export
#' ELEUTHIA_discriminative_report(
#'   glasso_fit = fit,
#'   cv_result = cv_result,
#'   selected = selected,
#'   final_model = final,
#'   coef_table = coefs,
#'   output_dir = "results/discriminative"
#' )
#'
#' }
ELEUTHIA_discriminative_report <- function(glasso_fit,
                                            cv_result,
                                            selected,
                                            final_model = NULL,
                                            coef_table = NULL,
                                            encoded = NULL,
                                            output_dir = NULL,
                                            prefix = "discriminative",
                                            save_tables = TRUE,
                                            save_plots = TRUE,
                                            plot_width = 10,
                                            plot_height = 8,
                                            verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!inherits(glasso_fit, "artemis_group_lasso")) {
    stop("glasso_fit must be an artemis_group_lasso object")
  }

  if (is.null(cv_result) || !is.list(cv_result)) {
    stop("cv_result must be the output from ARTEMIS_select_lambda()")
  }

  if (is.null(selected) || !is.list(selected)) {
    stop("selected must be the output from ARTEMIS_extract_selected_variables()")
  }

  # Create output directory if specified
  files_saved <- character(0)
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  }

  # Initialize report
  report <- character(0)
  tables <- list()

  # Helper functions
  add_line <- function(...) {
    report <<- c(report, paste0(...))
  }

  add_separator <- function(char = "=", width = 70) {
    add_line(paste(rep(char, width), collapse = ""))
  }

  # ---------------------------------------------------------------------------
  # Header
  # ---------------------------------------------------------------------------
  add_separator()
  add_line("DISCRIMINATIVE VARIABLE ANALYSIS REPORT")
  add_separator()
  add_line("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  add_line("")

  # ---------------------------------------------------------------------------
  # Section 1: Data Overview
  # ---------------------------------------------------------------------------
  add_separator("-", 70)
  add_line("1. DATA OVERVIEW")
  add_separator("-", 70)

  add_line("Samples: ", glasso_fit$n_samples)
  add_line("Features (encoded): ", glasso_fit$n_features)
  add_line("Variable groups: ", glasso_fit$n_groups)
  add_line("Response family: ", glasso_fit$family)
  add_line("Penalty type: ", glasso_fit$penalty)
  add_line("")

  if (!is.null(encoded)) {
    add_line("Encoding details:")
    add_line("  Target variable: ", encoded$target_name)
    n_quanti <- sum(encoded$var_mapping$var_types == "quantitative")
    n_cat <- sum(encoded$var_mapping$var_types == "categorical")
    n_ord <- sum(encoded$var_mapping$var_types == "ordinal")
    add_line("  Quantitative predictors: ", n_quanti)
    add_line("  Categorical predictors: ", n_cat)
    if (n_ord > 0) add_line("  Ordinal predictors: ", n_ord)
    add_line("  Encoding method: ", encoded$var_mapping$encoding)
    add_line("")
  }

  # ---------------------------------------------------------------------------
  # Section 2: Variable Selection
  # ---------------------------------------------------------------------------
  add_separator("-", 70)
  add_line("2. VARIABLE SELECTION (Group LASSO)")
  add_separator("-", 70)

  add_line("Lambda path range: [", round(min(glasso_fit$lambda), 4), ", ",
           round(max(glasso_fit$lambda), 4), "]")
  add_line("Number of lambda values: ", length(glasso_fit$lambda))
  add_line("")

  add_line("Cross-validation results (", cv_result$cv_fit$nfolds, "-fold):")
  add_line("  Lambda.min: ", round(cv_result$lambda.min, 4),
           " (", cv_result$n_selected_min, " groups selected)")
  add_line("  Lambda.1se: ", round(cv_result$lambda.1se, 4),
           " (", cv_result$n_selected_1se, " groups selected)")
  add_line("")

  add_line("Selection at lambda = ", round(selected$lambda, 4), ":")
  add_line("  Variables selected: ", selected$n_selected, " of ", selected$n_total)
  add_line("  Selection rate: ", round(100 * selected$n_selected / selected$n_total, 1), "%")
  add_line("")

  # List selected variables
  add_line("Selected variables:")
  for (i in seq_along(selected$selected_names)) {
    add_line("  ", i, ". ", selected$selected_names[i])
  }
  add_line("")

  # Create selected variables table
  selected_table <- data.frame(
    rank = seq_along(selected$selected_names),
    variable = selected$selected_names,
    group_id = selected$selected_groups,
    stringsAsFactors = FALSE
  )
  tables$selected_variables <- selected_table

  # ---------------------------------------------------------------------------
  # Section 3: Final Model
  # ---------------------------------------------------------------------------
  if (!is.null(final_model) && inherits(final_model, "artemis_final_model")) {
    add_separator("-", 70)
    add_line("3. FINAL MODEL (Unpenalized)")
    add_separator("-", 70)

    add_line("Family: ", final_model$family)
    add_line("Samples used: ", final_model$n_samples)
    add_line("Converged: ", final_model$converged)
    add_line("")

    # Model fit statistics
    model <- final_model$model
    if (final_model$family == "binomial") {
      add_line("Model fit:")
      add_line("  Null deviance: ", round(model$null.deviance, 2),
               " on ", model$df.null, " df")
      add_line("  Residual deviance: ", round(model$deviance, 2),
               " on ", model$df.residual, " df")
      add_line("  AIC: ", round(model$aic, 2))

      # Pseudo R-squared (McFadden)
      pseudo_r2 <- 1 - (model$deviance / model$null.deviance)
      add_line("  Pseudo R-squared (McFadden): ", round(pseudo_r2, 3))
    } else if (final_model$family == "multinomial") {
      add_line("Model fit:")
      add_line("  Residual deviance: ", round(model$deviance, 2))
      add_line("  AIC: ", round(model$AIC, 2))
      add_line("  Classes: ", paste(model$lev, collapse = ", "))
    } else if (final_model$family == "gaussian") {
      add_line("Model fit:")
      model_summary <- summary(model)
      residual_se <- sqrt(model_summary$dispersion)
      add_line("  Residual SE: ", round(residual_se, 4))
      r_sq <- 1 - (model$deviance / model$null.deviance)
      add_line("  R-squared: ", round(r_sq, 3))
    } else {
      # Poisson or other families
      add_line("Model fit:")
      add_line("  Residual deviance: ", round(model$deviance, 2), " on ",
               model$df.residual, " df")
      add_line("  AIC: ", round(model$aic, 2))
    }
    add_line("")
  }

  # ---------------------------------------------------------------------------
  # Section 4: Coefficients
  # ---------------------------------------------------------------------------
  if (!is.null(coef_table) && is.data.frame(coef_table)) {
    add_separator("-", 70)
    add_line("4. COEFFICIENTS")
    add_separator("-", 70)

    # Determine if odds ratios or coefficients
    is_or <- any(coef_table$estimate > 0) && min(coef_table$estimate, na.rm = TRUE) > 0
    est_label <- if (is_or) "Odds Ratio" else "Coefficient"

    # Check if multinomial (has 'class' column)
    is_multinomial <- "class" %in% colnames(coef_table)

    if (is_multinomial) {
      add_line("Multinomial model coefficients (", est_label, "s):")
      add_line("Reference class: ", final_model$model$lev[1])
      add_line("")

      for (cls in unique(coef_table$class)) {
        add_line("Class: ", cls, " vs reference")
        add_line(sprintf("  %-25s %10s %18s %10s",
                         "Variable", est_label, "95% CI", "P-value"))
        add_line(paste(rep("-", 65), collapse = ""))

        cls_data <- coef_table[coef_table$class == cls, ]
        for (i in seq_len(nrow(cls_data))) {
          row <- cls_data[i, ]
          ci_str <- sprintf("[%s, %s]", row$ci_lower, row$ci_upper)
          sig <- if (!is.null(row$significance)) row$significance else ""
          add_line(sprintf("  %-25s %10s %18s %10s %s",
                           substr(row$variable, 1, 25),
                           row$estimate,
                           ci_str,
                           row$p_value,
                           sig))
        }
        add_line("")
      }
    } else {
      add_line("Model coefficients (", est_label, "s):")
      add_line(sprintf("  %-25s %10s %18s %10s",
                       "Variable", est_label, "95% CI", "P-value"))
      add_line(paste(rep("-", 65), collapse = ""))

      for (i in seq_len(nrow(coef_table))) {
        row <- coef_table[i, ]
        ci_str <- sprintf("[%s, %s]", row$ci_lower, row$ci_upper)
        sig <- if (!is.null(row$significance)) row$significance else ""
        add_line(sprintf("  %-25s %10s %18s %10s %s",
                         substr(row$variable, 1, 25),
                         row$estimate,
                         ci_str,
                         row$p_value,
                         sig))
      }
      add_line("")
    }

    add_line("Signif. codes: *** p<0.001, ** p<0.01, * p<0.05")
    add_line("")

    tables$coefficients <- coef_table
  }

  # ---------------------------------------------------------------------------
  # Section 5: Interpretation Guide
  # ---------------------------------------------------------------------------
  add_separator("-", 70)
  add_line("5. INTERPRETATION GUIDE")
  add_separator("-", 70)

  if (glasso_fit$family == "binomial") {
    add_line("For binary outcomes (logistic regression):")
    add_line("  - Odds Ratio > 1: Higher predictor value -> higher probability of outcome")
    add_line("  - Odds Ratio < 1: Higher predictor value -> lower probability of outcome")
    add_line("  - Odds Ratio = 1: No association")
    add_line("  - If 95% CI includes 1: Not statistically significant at p=0.05")
  } else if (glasso_fit$family == "multinomial") {
    add_line("For multi-class outcomes (multinomial regression):")
    add_line("  - Each coefficient compares a class to the reference class")
    add_line("  - Odds Ratio > 1: Variable increases odds of that class vs reference")
    add_line("  - Interpret each class comparison separately")
  } else {
    add_line("For continuous outcomes (linear regression):")
    add_line("  - Coefficient: Change in outcome per unit change in predictor")
    add_line("  - Positive coefficient: Positive association")
    add_line("  - Negative coefficient: Negative association")
  }
  add_line("")

  add_line("Variable selection notes:")
  add_line("  - Variables were selected using Group LASSO with cross-validation")
  add_line("  - Lambda.1se provides more parsimonious selection (recommended)")
  add_line("  - Selected variables are consistently predictive across CV folds")
  add_line("  - Final model p-values should be interpreted cautiously (post-selection)")
  add_line("")

  add_separator()
  add_line("END OF REPORT")
  add_separator()

  # ---------------------------------------------------------------------------
  # Print to console
  # ---------------------------------------------------------------------------
  if (verbose) {
    cat("[ELEUTHIA] ",paste(report, collapse = "\n"), "\n")
  }

  # ---------------------------------------------------------------------------
  # Save files
  # ---------------------------------------------------------------------------
  if (!is.null(output_dir)) {

    # Save text report
    report_file <- file.path(output_dir, paste0(prefix, "_report.txt"))
    writeLines(report, report_file)
    files_saved <- c(files_saved, report_file)

    # Save tables
    if (save_tables) {
      for (table_name in names(tables)) {
        table_file <- file.path(output_dir, paste0(prefix, "_", table_name, ".csv"))
        write.csv(tables[[table_name]], table_file, row.names = FALSE)
        files_saved <- c(files_saved, table_file)
      }
    }

    # Save plots
    if (save_plots) {
      # Forest plot (only if coef_table available)
      if (!is.null(coef_table) && !("class" %in% colnames(coef_table))) {
        tryCatch({
          forest_plot <- AETHER_plot_coefficients(coef_table)
          forest_file <- file.path(output_dir, paste0(prefix, "_forest_plot.png"))
          ggplot2::ggsave(forest_file, forest_plot, width = plot_width, height = plot_height)
          files_saved <- c(files_saved, forest_file)
        }, error = function(e) {
          warning("Could not save forest plot: ", e$message)
        })
      }

      # Lambda path plot
      tryCatch({
        lambda_plot <- AETHER_plot_lambda_path(cv_result)
        lambda_file <- file.path(output_dir, paste0(prefix, "_lambda_path.png"))
        ggplot2::ggsave(lambda_file, lambda_plot, width = plot_width, height = plot_height)
        files_saved <- c(files_saved, lambda_file)
      }, error = function(e) {
        warning("Could not save lambda path plot: ", e$message)
      })
    }

    if (verbose) {
      cat("[ELEUTHIA] Files saved to:", output_dir, "\n")
      for (f in files_saved) {
        cat("[ELEUTHIA]   -", basename(f), "\n")
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------
  invisible(list(
    summary = report,
    tables = tables,
    files_saved = files_saved
  ))
}


#' Quick Summary of Discriminative Analysis
#'
#' @description Print a quick one-page summary without file export.
#'
#' @param glasso_fit An artemis_group_lasso object.
#' @param cv_result Output from ARTEMIS_select_lambda().
#' @param selected Output from ARTEMIS_extract_selected_variables().
#'
#' @return Invisibly returns NULL.
#'
#' @export
ELEUTHIA_discriminative_quick_summary <- function(glasso_fit,
                                                   cv_result,
                                                   selected) {

  cat("[ELEUTHIA] DISCRIMINATIVE ANALYSIS SUMMARY \n\n")

  cat("    DATA:\n")
  cat("        Samples:", glasso_fit$n_samples, "| Variables:", glasso_fit$n_groups, "\n")
  cat("        Family:", glasso_fit$family, "| Penalty:", glasso_fit$penalty, "\n\n")

  cat("    SELECTION:\n")
  cat("        Lambda.1se:", round(cv_result$lambda.1se, 4),
      "->", cv_result$n_selected_1se, "variables\n")
  cat("        Lambda.min:", round(cv_result$lambda.min, 4),
      "->", cv_result$n_selected_min, "variables\n\n")

  cat("    SELECTED VARIABLES:\n")
  n_show <- min(10, length(selected$selected_names))
  for (i in seq_len(n_show)) {
    cat("        ", i, ". ", selected$selected_names[i], "\n", sep = "")
  }
  if (length(selected$selected_names) > 10) {
    cat("        ... and", length(selected$selected_names) - 10, "more\n")
  }
  cat("\n")

  cat("=============================================\n")

  invisible(NULL)
}
