#' Quality Control Assessment for Mixed Data
#'
#' @description Assess a mixed (quantitative + qualitative) dataset and report
#' QC metrics to inform filtering decisions. This function reports and flags
#' potential issues but does NOT automatically remove anything - the user
#' reviews the report and decides what to filter.
#'
#' @param data Data frame with mixed variable types.
#' @param na_threshold Numeric (0-1). Variables/samples with NA proportion
#'   exceeding this threshold will be flagged. Default = 0.2 (20%).
#' @param variance_threshold Numeric. Quantitative variables with variance
#'   below this threshold will be flagged as near-zero variance. Default = 1e-10.
#' @param dominance_threshold Numeric (0-1). Qualitative variables where one
#'   level represents more than this proportion will be flagged as near-constant.
#'   Default = 0.99 (99%).
#' @param correlation_threshold Numeric (0-1). Pairs of quantitative variables
#'   with absolute correlation above this threshold will be reported.
#'   Default = 0.9.
#' @param check_correlations Logical. Whether to compute correlation matrix
#'   for quantitative variables. Default = TRUE.
#' @param check_samples Logical. Whether to compute sample-level QC metrics.
#'   Default = TRUE.
#' @param check_type_issues Logical. Whether to check for type inconsistencies
#'   and NA-like strings. Default = TRUE.
#' @param na_strings Character vector. Strings to detect as potential NA values.
#'   Default includes common patterns: "N/A", "NA", "None", "NULL", "", etc.
#'   Set to NULL to disable NA-string detection.
#' @param numeric_threshold Numeric (0-1). For character columns, if this
#'   proportion of non-NA values can be parsed as numeric, flag as potential
#'   type issue. Default = 0.5 (50%).
#' @param categorical_threshold Integer. Numeric columns with this many or fewer
#'   unique values will be flagged as potentially categorical (might need to be
#'   converted to factor for FAMD). Default = 10.
#' @param sample_id_col Character or integer. Column name or index containing
#'   sample IDs. Recommended to specify this for meaningful sample identification.
#'   When specified, the column is validated for uniqueness, missing values, and
#'   empty strings. If NULL (default), row names or row numbers are used.
#' @param verbose Logical. Print progress messages. Default = TRUE.
#'
#' @return A list with class "hades_qc" containing:
#' \describe{
#'   \item{variables}{Data frame with per-variable QC metrics and flags}
#'   \item{samples}{Data frame with per-sample QC metrics and flags (if check_samples = TRUE)}
#'   \item{correlations}{Data frame of highly correlated variable pairs (if check_correlations = TRUE)}
#'   \item{type_issues}{List with type inconsistencies, NA-like strings, and
#'     special characters found (if check_type_issues = TRUE). Includes:
#'     potential_numeric, potential_categorical, na_like_strings,
#'     mixed_na_representations, special_chars_in_values}
#'   \item{id_issues}{List with sample ID validation results (if sample_id_col specified)}
#'   \item{summary}{List with counts of flagged variables/samples}
#'   \item{thresholds}{List of thresholds used for flagging}
#' }
#'
#' @details
#' This function is designed for use before FAMD or other mixed-data analyses.
#' It helps identify:
#' \itemize{
#'   \item Variables with excessive missing data
#'   \item Variables with near-zero variance (uninformative)
#'   \item Highly correlated variable pairs (redundancy)
#'   \item Samples with excessive missing data
#'   \item Type inconsistencies (character columns that look numeric)
#'   \item NA-like strings that weren't parsed as NA (e.g., "N/A", "None")
#'   \item Numeric columns that may be categorical (few unique values)
#'   \item Special characters in categorical values (spaces, underscores) that
#'         may cause issues with FAMD category naming
#' }
#'
#' The user should review the output and make informed decisions about what
#' to filter before proceeding with analysis.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic QC assessment
#' qc_report <- HADES_qc_mixed_data(my_data)
#' print(qc_report)
#'
#' # View flagged variables
#' qc_report$variables[qc_report$variables$flagged_any, ]
#'
#' # View highly correlated pairs
#' qc_report$correlations
#'
#' # View type issues (NA-like strings, potential numeric columns)
#' qc_report$type_issues
#'
#' }
HADES_qc_mixed_data <- function(data,
                                 na_threshold = 0.2,
                                 variance_threshold = 1e-10,
                                 dominance_threshold = 0.99,
                                 correlation_threshold = 0.9,
                                 check_correlations = TRUE,
                                 check_samples = TRUE,
                                 check_type_issues = TRUE,
                                 na_strings = c(
                                   # Standard NA representations
                                   "NA", "N/A", "na", "n/a", "Na", "N/a",
                                   # None variants
                                   "None", "none", "NONE",
                                   # NULL variants
                                   "NULL", "null", "Null",
                                   # NaN variants
                                   "NaN", "nan", "NAN",
                                   # Empty/placeholder
                                   "", " ", "  ", ".", "-", "--", "---",
                                   # Missing variants
                                   "missing", "Missing", "MISSING",
                                   # Not available variants
                                   "not available", "Not Available", "NOT AVAILABLE",
                                   # Unknown variants
                                   "unknown", "Unknown", "UNKNOWN", "unk", "UNK",
                                   # Excel errors
                                   "#N/A", "#NA", "#VALUE!", "#REF!", "#DIV/0!", "#NAME?",
                                   # Other common patterns
                                   "n.a.", "N.A.", "n.a", "N.A",
                                   "nil", "NIL", "Nil",
                                   "undefined", "UNDEFINED"
                                 ),
                                 numeric_threshold = 0.5,
                                 categorical_threshold = 10,
                                 sample_id_col = NULL,
                                 verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  if (nrow(data) == 0) {
    stop("Data has no rows")
  }
  if (ncol(data) == 0) {
    stop("Data has no columns")
  }

  # Validate thresholds
  if (na_threshold < 0 || na_threshold > 1) {
    stop("na_threshold must be between 0 and 1")
  }
  if (dominance_threshold < 0 || dominance_threshold > 1) {
    stop("dominance_threshold must be between 0 and 1")
  }
  if (correlation_threshold < 0 || correlation_threshold > 1) {
    stop("correlation_threshold must be between 0 and 1")
  }

  if (verbose) {
    cat("[HADES] Mixed Data QC Assessment\n")
    cat("===============================\n")
    cat("    Input:", nrow(data), "samples,", ncol(data), "variables\n\n")
  }

  # ---------------------------------------------------------------------------
  # Step 1: Identify variable types
  # ---------------------------------------------------------------------------
  var_types <- sapply(data, function(x) {
    if (is.numeric(x)) "quantitative"
    else if (is.factor(x) || is.character(x) || is.logical(x)) "qualitative"
    else "other"
  })

  n_quanti <- sum(var_types == "quantitative")
  n_quali <- sum(var_types == "qualitative")
  n_other <- sum(var_types == "other")

  if (verbose) {
    cat("[HADES] Variable types detected:\n")
    cat("    Quantitative:", n_quanti, "\n")
    cat("    Qualitative:", n_quali, "\n")
    if (n_other > 0) cat("    Other:", n_other, "(excluded from some checks)\n")
    cat("\n")
  }

  # ---------------------------------------------------------------------------
  # Step 2: Variable-level QC
  # ---------------------------------------------------------------------------
  if (verbose) cat("[HADES] Assessing variables...\n")

  var_report <- data.frame(
    variable = names(data),
    type = var_types,
    stringsAsFactors = FALSE
  )

  # Missing data
  var_report$n_total <- nrow(data)
  var_report$n_missing <- sapply(data, function(x) sum(is.na(x)))
  var_report$pct_missing <- var_report$n_missing / var_report$n_total

  # Variance / level info
  var_report$variance <- NA_real_
  var_report$n_levels <- NA_integer_
  var_report$max_level_pct <- NA_real_

  for (i in seq_along(data)) {
    col <- data[[i]]
    col_complete <- col[!is.na(col)]

    if (var_types[i] == "quantitative") {
      if (length(col_complete) > 1) {
        var_report$variance[i] <- var(col_complete)
      } else {
        var_report$variance[i] <- 0
      }
    } else if (var_types[i] == "qualitative") {
      if (length(col_complete) > 0) {
        level_table <- table(col_complete)
        var_report$n_levels[i] <- length(level_table)
        var_report$max_level_pct[i] <- max(level_table) / sum(level_table)
      } else {
        var_report$n_levels[i] <- 0
        var_report$max_level_pct[i] <- NA
      }
    }
  }

  # Flagging
  var_report$flagged_na <- var_report$pct_missing > na_threshold

  var_report$flagged_variance <- FALSE
  var_report$flagged_variance[var_types == "quantitative"] <-
    var_report$variance[var_types == "quantitative"] < variance_threshold

  var_report$flagged_dominance <- FALSE
  var_report$flagged_dominance[var_types == "qualitative"] <-
    var_report$n_levels[var_types == "qualitative"] <= 1 |
    var_report$max_level_pct[var_types == "qualitative"] > dominance_threshold

  # Handle NAs in flags (from completely missing columns)
  var_report$flagged_variance[is.na(var_report$flagged_variance)] <- TRUE
  var_report$flagged_dominance[is.na(var_report$flagged_dominance)] <- TRUE

  var_report$flagged_any <- var_report$flagged_na |
    var_report$flagged_variance |
    var_report$flagged_dominance

  # ---------------------------------------------------------------------------
  # Step 3: Type consistency and NA-like string detection
  # ---------------------------------------------------------------------------
  type_issues <- NULL

  if (check_type_issues) {
    if (verbose) cat("[HADES] Checking for type inconsistencies and NA-like strings...\n")

    # Initialize results
    potential_numeric <- list()
    potential_categorical <- list()
    na_like_found <- list()
    mixed_na_representations <- list()
    special_chars_in_values <- list()  # Spaces and special chars in factor levels

    for (i in seq_along(data)) {
      col <- data[[i]]
      col_name <- names(data)[i]

      # ----- Check numeric columns for potential categorical -----
      if (is.numeric(col)) {
        col_complete <- col[!is.na(col)]
        if (length(col_complete) > 0) {
          unique_vals <- unique(col_complete)
          n_unique <- length(unique_vals)

          if (n_unique <= categorical_threshold) {
            # Check if values look like codes (integers or simple decimals)
            all_integer_like <- all(col_complete == floor(col_complete))

            potential_categorical[[col_name]] <- list(
              variable = col_name,
              n_unique = n_unique,
              unique_values = sort(unique_vals),
              all_integer = all_integer_like
            )
          }
        }
      }

      # Only check character/factor columns for NA-like strings and potential numeric
      if (is.character(col) || is.factor(col)) {
        col_char <- as.character(col)
        col_non_na <- col_char[!is.na(col_char)]

        if (length(col_non_na) == 0) next

        # ----- Check for special characters in values -----
        # These can cause issues with FAMD category naming (e.g., "VarName_Category")
        unique_vals <- unique(col_non_na)
        # Check for spaces
        has_space <- grepl(" ", unique_vals, fixed = TRUE)
        # Check for other problematic characters (underscore is used in FAMD naming)
        has_underscore <- grepl("_", unique_vals, fixed = TRUE)
        # Check for other special chars that might cause issues
        has_special <- grepl("[^a-zA-Z0-9._-]", unique_vals) & !has_space

        if (any(has_space) || any(has_underscore) || any(has_special)) {
          problem_vals <- list()
          if (any(has_space)) {
            problem_vals$spaces <- unique_vals[has_space]
          }
          if (any(has_underscore)) {
            problem_vals$underscores <- unique_vals[has_underscore]
          }
          if (any(has_special)) {
            problem_vals$special <- unique_vals[has_special]
          }
          special_chars_in_values[[col_name]] <- problem_vals
        }

        # ----- Check for NA-like strings -----
        if (!is.null(na_strings) && length(na_strings) > 0) {
          # Trim whitespace for comparison
          col_trimmed <- trimws(col_non_na)

          # Find matches (case-sensitive first, then we report what was found)
          na_matches <- col_trimmed %in% na_strings

          if (any(na_matches)) {
            matched_values <- unique(col_trimmed[na_matches])
            matched_counts <- table(col_trimmed[na_matches])

            na_like_found[[col_name]] <- data.frame(
              variable = col_name,
              na_string = names(matched_counts),
              count = as.integer(matched_counts),
              stringsAsFactors = FALSE
            )

            # Check for multiple different NA representations in same column
            if (length(matched_values) > 1) {
              mixed_na_representations[[col_name]] <- matched_values
            }
          }
        }

        # ----- Check if column looks numeric -----
        # Remove NA-like strings before checking
        col_for_numeric_check <- col_non_na
        if (!is.null(na_strings)) {
          col_for_numeric_check <- col_non_na[!trimws(col_non_na) %in% na_strings]
        }

        if (length(col_for_numeric_check) > 0) {
          # Try to parse as numeric
          parsed <- suppressWarnings(as.numeric(col_for_numeric_check))
          n_parseable <- sum(!is.na(parsed))
          pct_parseable <- n_parseable / length(col_for_numeric_check)

          if (pct_parseable >= numeric_threshold && n_parseable > 0) {
            # Find the values that prevented numeric parsing
            non_numeric_values <- unique(col_for_numeric_check[is.na(parsed)])

            potential_numeric[[col_name]] <- list(
              variable = col_name,
              pct_numeric = pct_parseable,
              n_numeric = n_parseable,
              n_non_numeric = length(col_for_numeric_check) - n_parseable,
              non_numeric_values = non_numeric_values
            )
          }
        }
      }
    }

    # Build type_issues output
    type_issues <- list(
      potential_numeric = potential_numeric,
      potential_categorical = potential_categorical,
      na_like_strings = na_like_found,
      mixed_na_representations = mixed_na_representations,
      special_chars_in_values = special_chars_in_values
    )

    # Add flags to var_report
    var_report$has_na_strings <- names(data) %in% names(na_like_found)
    var_report$potential_numeric <- names(data) %in% names(potential_numeric)
    var_report$potential_categorical <- names(data) %in% names(potential_categorical)
    var_report$has_special_chars <- names(data) %in% names(special_chars_in_values)
    var_report$flagged_type_issue <- var_report$has_na_strings |
      var_report$potential_numeric |
      var_report$potential_categorical |
      var_report$has_special_chars

    # Update flagged_any
    var_report$flagged_any <- var_report$flagged_any | var_report$flagged_type_issue
  } else {
    # Add empty columns if not checking
    var_report$has_na_strings <- FALSE
    var_report$potential_numeric <- FALSE
    var_report$potential_categorical <- FALSE
    var_report$has_special_chars <- FALSE
    var_report$flagged_type_issue <- FALSE
  }

  # ---------------------------------------------------------------------------
  # Step 4: Correlation analysis (quantitative variables only)
  # ---------------------------------------------------------------------------
  cor_report <- NULL

  if (check_correlations && n_quanti >= 2) {
    if (verbose) cat("[HADES] Computing correlations for quantitative variables...\n")

    quanti_cols <- names(var_types)[var_types == "quantitative"]
    quanti_data <- data[, quanti_cols, drop = FALSE]

    # Compute correlation matrix (pairwise complete)
    cor_matrix <- cor(quanti_data, use = "pairwise.complete.obs")

    # Find pairs above threshold
    cor_pairs <- list()
    for (i in 1:(ncol(cor_matrix) - 1)) {
      for (j in (i + 1):ncol(cor_matrix)) {
        r <- cor_matrix[i, j]
        if (!is.na(r) && abs(r) > correlation_threshold) {
          cor_pairs[[length(cor_pairs) + 1]] <- data.frame(
            var1 = colnames(cor_matrix)[i],
            var2 = colnames(cor_matrix)[j],
            correlation = r,
            stringsAsFactors = FALSE
          )
        }
      }
    }

    if (length(cor_pairs) > 0) {
      cor_report <- do.call(rbind, cor_pairs)
      cor_report <- cor_report[order(-abs(cor_report$correlation)), ]
      rownames(cor_report) <- NULL
    } else {
      cor_report <- data.frame(
        var1 = character(0),
        var2 = character(0),
        correlation = numeric(0),
        stringsAsFactors = FALSE
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Step 5: Sample-level QC
  # ---------------------------------------------------------------------------
  sample_report <- NULL
  id_issues <- NULL

  if (check_samples) {
    if (verbose) cat("[HADES] Assessing samples...\n")

    # Determine sample IDs and validate if column specified
    id_source <- "generated"  # Track where IDs came from

    if (!is.null(sample_id_col)) {
      # Validate that the column exists
      if (is.character(sample_id_col)) {
        if (!sample_id_col %in% names(data)) {
          stop("sample_id_col '", sample_id_col, "' not found in data")
        }
        sample_ids <- data[[sample_id_col]]
        id_col_name <- sample_id_col
      } else {
        if (sample_id_col < 1 || sample_id_col > ncol(data)) {
          stop("sample_id_col index ", sample_id_col, " is out of range")
        }
        sample_ids <- data[[sample_id_col]]
        id_col_name <- names(data)[sample_id_col]
      }
      id_source <- "column"

      # ----- Validate ID column -----
      id_issues <- list(
        column = id_col_name,
        n_missing = sum(is.na(sample_ids)),
        missing_rows = which(is.na(sample_ids)),
        n_duplicates = 0,
        duplicate_values = character(0),
        duplicate_rows = list(),
        n_empty_strings = 0,
        empty_string_rows = integer(0)
      )

      # Check for missing IDs
      if (id_issues$n_missing > 0 && verbose) {
        cat("[HADES] WARNING:", id_issues$n_missing, "samples have missing IDs\n")
      }

      # Check for empty strings (common issue)
      sample_ids_char <- as.character(sample_ids)
      empty_mask <- !is.na(sample_ids_char) & trimws(sample_ids_char) == ""
      id_issues$n_empty_strings <- sum(empty_mask)
      id_issues$empty_string_rows <- which(empty_mask)

      if (id_issues$n_empty_strings > 0 && verbose) {
        cat("[HADES] WARNING:", id_issues$n_empty_strings, "samples have empty string IDs\n")
      }

      # Check for duplicates (only among non-NA, non-empty values)
      valid_ids <- sample_ids_char[!is.na(sample_ids_char) & trimws(sample_ids_char) != ""]
      dup_table <- table(valid_ids)
      duplicated_vals <- names(dup_table)[dup_table > 1]

      if (length(duplicated_vals) > 0) {
        id_issues$n_duplicates <- length(duplicated_vals)
        id_issues$duplicate_values <- duplicated_vals

        # Find which rows have each duplicate
        for (dup_val in duplicated_vals) {
          id_issues$duplicate_rows[[dup_val]] <- which(sample_ids_char == dup_val)
        }

        if (verbose) {
          cat("[HADES] WARNING:", id_issues$n_duplicates, "duplicate ID value(s) found\n")
        }
      }

      # Convert sample_ids to character for consistency
      sample_ids <- sample_ids_char

    } else if (!is.null(rownames(data)) && !identical(rownames(data), as.character(1:nrow(data)))) {
      sample_ids <- rownames(data)
      id_source <- "rownames"
    } else {
      sample_ids <- paste0("row_", seq_len(nrow(data)))
      id_source <- "generated"
      if (verbose) {
        cat("  NOTE: No sample_id_col specified. Using row numbers as IDs.\n")
        cat("        Specify sample_id_col for meaningful sample identification.\n")
      }
    }

    sample_report <- data.frame(
      sample_id = sample_ids,
      n_total = ncol(data),
      id_source = id_source,
      stringsAsFactors = FALSE
    )

    # Exclude ID column from missing count if it was specified
    if (!is.null(sample_id_col)) {
      data_for_missing <- data[, names(data) != id_col_name, drop = FALSE]
      sample_report$n_total <- ncol(data_for_missing)
      sample_report$n_missing <- rowSums(is.na(data_for_missing))
    } else {
      sample_report$n_missing <- rowSums(is.na(data))
    }
    sample_report$pct_missing <- sample_report$n_missing / sample_report$n_total
    sample_report$flagged_na <- sample_report$pct_missing > na_threshold
  }

  # ---------------------------------------------------------------------------
  # Step 6: Summary
  # ---------------------------------------------------------------------------
  summary_stats <- list(
    n_variables = ncol(data),
    n_samples = nrow(data),
    n_quanti = n_quanti,
    n_quali = n_quali,
    variables_flagged_na = sum(var_report$flagged_na),
    variables_flagged_variance = sum(var_report$flagged_variance),
    variables_flagged_dominance = sum(var_report$flagged_dominance),
    variables_flagged_any = sum(var_report$flagged_any)
  )

  if (check_samples) {
    summary_stats$samples_flagged_na <- sum(sample_report$flagged_na)
    summary_stats$id_source <- sample_report$id_source[1]

    if (!is.null(id_issues)) {
      summary_stats$id_missing <- id_issues$n_missing
      summary_stats$id_empty_strings <- id_issues$n_empty_strings
      summary_stats$id_duplicates <- id_issues$n_duplicates
      summary_stats$id_issues_any <- (id_issues$n_missing + id_issues$n_empty_strings + id_issues$n_duplicates) > 0
    }
  }

  if (check_correlations && !is.null(cor_report)) {
    summary_stats$n_high_correlation_pairs <- nrow(cor_report)
  }

  if (check_type_issues && !is.null(type_issues)) {
    summary_stats$variables_with_na_strings <- sum(var_report$has_na_strings)
    summary_stats$variables_potential_numeric <- sum(var_report$potential_numeric)
    summary_stats$variables_potential_categorical <- sum(var_report$potential_categorical)
    summary_stats$n_mixed_na_representations <- length(type_issues$mixed_na_representations)
    summary_stats$variables_with_special_chars <- sum(var_report$has_special_chars)
  }

  # ---------------------------------------------------------------------------
  # Build result
  # ---------------------------------------------------------------------------
  result <- list(
    variables = var_report,
    samples = sample_report,
    correlations = cor_report,
    type_issues = type_issues,
    id_issues = id_issues,
    summary = summary_stats,
    thresholds = list(
      na_threshold = na_threshold,
      variance_threshold = variance_threshold,
      dominance_threshold = dominance_threshold,
      correlation_threshold = correlation_threshold,
      numeric_threshold = numeric_threshold,
      categorical_threshold = categorical_threshold,
      na_strings = na_strings
    )
  )

  class(result) <- c("hades_qc", "list")

  if (verbose) {
    cat("\n")
    cat("[HADES] QC Summary:\n")
    cat("-----------\n")
    cat("    Variables flagged (high NA):", summary_stats$variables_flagged_na, "\n")
    cat("    Variables flagged (low variance):", summary_stats$variables_flagged_variance, "\n")
    cat("    Variables flagged (dominant level):", summary_stats$variables_flagged_dominance, "\n")
    if (check_type_issues) {
      cat("    Variables with NA-like strings:", summary_stats$variables_with_na_strings, "\n")
      cat("    Variables potentially numeric:", summary_stats$variables_potential_numeric, "\n")
      cat("    Variables potentially categorical:", summary_stats$variables_potential_categorical, "\n")
      cat("    Variables with special chars in values:", summary_stats$variables_with_special_chars, "\n")
    }
    cat("    Variables flagged (any reason):", summary_stats$variables_flagged_any, "\n")
    if (check_samples) {
      cat("    Samples flagged (high NA):", summary_stats$samples_flagged_na, "\n")
      if (!is.null(id_issues) && summary_stats$id_issues_any) {
        cat("    Sample ID issues found: ")
        issues <- c()
        if (id_issues$n_missing > 0) issues <- c(issues, paste0(id_issues$n_missing, " missing"))
        if (id_issues$n_empty_strings > 0) issues <- c(issues, paste0(id_issues$n_empty_strings, " empty"))
        if (id_issues$n_duplicates > 0) issues <- c(issues, paste0(id_issues$n_duplicates, " duplicate values"))
        cat(paste(issues, collapse = ", "), "\n")
      }
    }
    if (check_correlations && !is.null(cor_report)) {
      cat("    High correlation pairs:", summary_stats$n_high_correlation_pairs, "\n")
    }
    cat("\n    Use print() for detailed report or access $variables, $samples, $correlations, $type_issues, $id_issues directly.\n")
  }

  return(result)
}


#' Print method for HADES QC results
#'
#' @param x A hades_qc object
#' @param show_flagged_only Logical. Only show flagged items. Default = TRUE.
#' @param ... Additional arguments (ignored)
#'
#' @method print hades_qc
#' @export
print.hades_qc <- function(x, show_flagged_only = TRUE, ...) {
  cat("HADES Mixed Data QC Report\n")
  cat("==========================\n\n")

  cat("Thresholds used:\n")
  cat("  NA threshold:", x$thresholds$na_threshold * 100, "%\n")
  cat("  Variance threshold:", x$thresholds$variance_threshold, "\n")
  cat("  Dominance threshold:", x$thresholds$dominance_threshold * 100, "%\n")
  cat("  Correlation threshold:", x$thresholds$correlation_threshold, "\n")
  cat("\n")

  cat("Summary:\n")
  cat("  Total variables:", x$summary$n_variables,
      "(", x$summary$n_quanti, "quanti,", x$summary$n_quali, "quali)\n")
  cat("  Total samples:", x$summary$n_samples, "\n")
  cat("\n")

  # Variables
  cat("VARIABLE FLAGS:\n")
  cat("---------------\n")
  flagged_vars <- x$variables[x$variables$flagged_any, ]


  if (nrow(flagged_vars) == 0) {
    cat("  No variables flagged.\n")
  } else {
    cat("  ", nrow(flagged_vars), "of", x$summary$n_variables, "variables flagged:\n\n")

    # High NA
    high_na <- x$variables[x$variables$flagged_na, ]
    if (nrow(high_na) > 0) {
      cat("  High missing data (>", x$thresholds$na_threshold * 100, "%):\n", sep = "")
      for (i in seq_len(nrow(high_na))) {
        cat("    -", high_na$variable[i], ":",
            round(high_na$pct_missing[i] * 100, 1), "% NA\n")
      }
      cat("\n")
    }

    # Low variance (quanti)
    low_var <- x$variables[x$variables$flagged_variance & x$variables$type == "quantitative", ]
    if (nrow(low_var) > 0) {
      cat("  Near-zero variance (quantitative):\n")
      for (i in seq_len(nrow(low_var))) {
        cat("    -", low_var$variable[i], ": var =",
            format(low_var$variance[i], scientific = TRUE, digits = 2), "\n")
      }
      cat("\n")
    }

    # Dominant level (quali)
    dom_level <- x$variables[x$variables$flagged_dominance & x$variables$type == "qualitative", ]
    if (nrow(dom_level) > 0) {
      cat("  Near-constant (qualitative):\n")
      for (i in seq_len(nrow(dom_level))) {
        if (dom_level$n_levels[i] <= 1) {
          cat("    -", dom_level$variable[i], ": only", dom_level$n_levels[i], "level(s)\n")
        } else {
          cat("    -", dom_level$variable[i], ":",
              round(dom_level$max_level_pct[i] * 100, 1), "% one level\n")
        }
      }
      cat("\n")
    }
  }

  # Samples
  if (!is.null(x$samples)) {
    cat("SAMPLE FLAGS:\n")
    cat("-------------\n")

    # Report ID source
    id_source <- x$samples$id_source[1]
    if (id_source == "generated") {
      cat("  NOTE: Sample IDs are row numbers (no sample_id_col specified)\n\n")
    } else if (id_source == "column") {
      cat("  Sample IDs from column:", x$id_issues$column, "\n\n")
    }

    flagged_samples <- x$samples[x$samples$flagged_na, ]

    if (nrow(flagged_samples) == 0) {
      cat("  No samples flagged for missing data.\n")
    } else {
      cat("  ", nrow(flagged_samples), "of", x$summary$n_samples, "samples flagged:\n\n")
      cat("  High missing data (>", x$thresholds$na_threshold * 100, "%):\n", sep = "")

      # Show first 10 if many
      n_show <- min(10, nrow(flagged_samples))
      for (i in seq_len(n_show)) {
        cat("    -", flagged_samples$sample_id[i], ":",
            round(flagged_samples$pct_missing[i] * 100, 1), "% NA\n")
      }
      if (nrow(flagged_samples) > 10) {
        cat("    ... and", nrow(flagged_samples) - 10, "more\n")
      }
    }
    cat("\n")
  }

  # ID Issues
  if (!is.null(x$id_issues)) {
    has_issues <- (x$id_issues$n_missing + x$id_issues$n_empty_strings + x$id_issues$n_duplicates) > 0

    if (has_issues) {
      cat("SAMPLE ID ISSUES:\n")
      cat("-----------------\n")
      cat("  ID column:", x$id_issues$column, "\n\n")

      # Missing IDs
      if (x$id_issues$n_missing > 0) {
        cat("  Missing IDs:", x$id_issues$n_missing, "sample(s)\n")
        n_show <- min(10, length(x$id_issues$missing_rows))
        cat("    Rows:", paste(head(x$id_issues$missing_rows, n_show), collapse = ", "))
        if (length(x$id_issues$missing_rows) > 10) {
          cat(" ... and", length(x$id_issues$missing_rows) - 10, "more")
        }
        cat("\n\n")
      }

      # Empty string IDs
      if (x$id_issues$n_empty_strings > 0) {
        cat("  Empty string IDs:", x$id_issues$n_empty_strings, "sample(s)\n")
        n_show <- min(10, length(x$id_issues$empty_string_rows))
        cat("    Rows:", paste(head(x$id_issues$empty_string_rows, n_show), collapse = ", "))
        if (length(x$id_issues$empty_string_rows) > 10) {
          cat(" ... and", length(x$id_issues$empty_string_rows) - 10, "more")
        }
        cat("\n\n")
      }

      # Duplicate IDs
      if (x$id_issues$n_duplicates > 0) {
        cat("  Duplicate IDs:", x$id_issues$n_duplicates, "value(s) appear multiple times\n\n")
        n_shown <- 0
        for (dup_val in names(x$id_issues$duplicate_rows)) {
          if (n_shown >= 10) {
            remaining <- length(x$id_issues$duplicate_rows) - n_shown
            cat("    ... and", remaining, "more duplicate values\n")
            break
          }
          rows <- x$id_issues$duplicate_rows[[dup_val]]
          cat("    -", paste0("\"", dup_val, "\""), "appears in rows:",
              paste(rows, collapse = ", "), "\n")
          n_shown <- n_shown + 1
        }
        cat("\n")
      }
    }
  }

  # Correlations
  if (!is.null(x$correlations) && nrow(x$correlations) > 0) {
    cat("HIGH CORRELATIONS:\n")
    cat("------------------\n")
    cat("  ", nrow(x$correlations), "pairs with |r| >", x$thresholds$correlation_threshold, ":\n\n")

    n_show <- min(10, nrow(x$correlations))
    for (i in seq_len(n_show)) {
      cat("    -", x$correlations$var1[i], "<->", x$correlations$var2[i],
          ": r =", round(x$correlations$correlation[i], 3), "\n")
    }
    if (nrow(x$correlations) > 10) {
      cat("    ... and", nrow(x$correlations) - 10, "more\n")
    }
    cat("\n")
  }

  # Type Issues
  if (!is.null(x$type_issues)) {
    has_na_strings <- length(x$type_issues$na_like_strings) > 0
    has_numeric <- length(x$type_issues$potential_numeric) > 0
    has_categorical <- length(x$type_issues$potential_categorical) > 0
    has_mixed <- length(x$type_issues$mixed_na_representations) > 0
    has_special <- length(x$type_issues$special_chars_in_values) > 0

    if (has_na_strings || has_numeric || has_categorical || has_mixed || has_special) {
      cat("TYPE ISSUES:\n")
      cat("------------\n")

      # NA-like strings found
      if (has_na_strings) {
        cat("  NA-like strings detected (not parsed as NA):\n\n")
        n_shown <- 0
        for (var_name in names(x$type_issues$na_like_strings)) {
          if (n_shown >= 10) {
            remaining <- length(x$type_issues$na_like_strings) - n_shown
            cat("    ... and", remaining, "more variables\n")
            break
          }
          df <- x$type_issues$na_like_strings[[var_name]]
          values_str <- paste(
            sapply(seq_len(nrow(df)), function(i) {
              val <- df$na_string[i]
              if (val == "") val <- '""'
              else if (trimws(val) == "") val <- paste0('"', val, '" (whitespace)')
              paste0(val, " (n=", df$count[i], ")")
            }),
            collapse = ", "
          )
          cat("    -", var_name, ":", values_str, "\n")
          n_shown <- n_shown + 1
        }
        cat("\n")
      }

      # Mixed NA representations (multiple different NA-like values in same column)
      if (has_mixed) {
        cat("  Columns with MULTIPLE NA representations (inconsistent):\n\n")
        for (var_name in names(x$type_issues$mixed_na_representations)) {
          values <- x$type_issues$mixed_na_representations[[var_name]]
          values_display <- sapply(values, function(v) {
            if (v == "") '""'
            else if (trimws(v) == "") paste0('"', v, '"')
            else v
          })
          cat("    -", var_name, ":", paste(values_display, collapse = ", "), "\n")
        }
        cat("\n")
      }

      # Potential numeric columns
      if (has_numeric) {
        cat("  Character columns that may be numeric:\n\n")
        n_shown <- 0
        for (var_name in names(x$type_issues$potential_numeric)) {
          if (n_shown >= 10) {
            remaining <- length(x$type_issues$potential_numeric) - n_shown
            cat("    ... and", remaining, "more variables\n")
            break
          }
          info <- x$type_issues$potential_numeric[[var_name]]
          cat("    -", var_name, ":",
              round(info$pct_numeric * 100, 1), "% parseable as numeric\n")
          if (length(info$non_numeric_values) > 0) {
            non_num_display <- head(info$non_numeric_values, 5)
            non_num_display <- sapply(non_num_display, function(v) {
              if (nchar(v) > 20) paste0(substr(v, 1, 17), "...")
              else v
            })
            cat("      Non-numeric values:", paste(non_num_display, collapse = ", "))
            if (length(info$non_numeric_values) > 5) {
              cat(" ... and", length(info$non_numeric_values) - 5, "more")
            }
            cat("\n")
          }
          n_shown <- n_shown + 1
        }
        cat("\n")
      }

      # Potential categorical columns (numeric with few unique values)
      if (has_categorical) {
        cat("  Numeric columns that may be categorical (<=",
            x$thresholds$categorical_threshold, "unique values):\n\n", sep = "")
        n_shown <- 0
        for (var_name in names(x$type_issues$potential_categorical)) {
          if (n_shown >= 15) {
            remaining <- length(x$type_issues$potential_categorical) - n_shown
            cat("    ... and", remaining, "more variables\n")
            break
          }
          info <- x$type_issues$potential_categorical[[var_name]]
          # Format unique values for display
          vals_display <- info$unique_values
          if (length(vals_display) > 8) {
            vals_display <- c(head(vals_display, 4), "...", tail(vals_display, 3))
          }
          vals_str <- paste(vals_display, collapse = ", ")

          cat("    -", var_name, ":", info$n_unique, "unique values")
          if (info$all_integer) cat(" [integers]")
          cat("\n")
          cat("      Values:", vals_str, "\n")
          n_shown <- n_shown + 1
        }
        cat("\n")
        cat("  TIP: Consider converting these to factors with as.factor() before FAMD\n")
        cat("       if they represent categories rather than continuous measurements.\n\n")
      }

      # Special characters in values (spaces, underscores, etc.)
      if (has_special) {
        cat("  Values with special characters (may cause FAMD issues):\n\n")
        n_shown <- 0
        for (var_name in names(x$type_issues$special_chars_in_values)) {
          if (n_shown >= 15) {
            remaining <- length(x$type_issues$special_chars_in_values) - n_shown
            cat("    ... and", remaining, "more variables\n")
            break
          }
          issues <- x$type_issues$special_chars_in_values[[var_name]]

          # Report spaces (most common issue)
          if (!is.null(issues$spaces)) {
            space_examples <- head(issues$spaces, 3)
            space_str <- paste0("\"", space_examples, "\"", collapse = ", ")
            if (length(issues$spaces) > 3) {
              space_str <- paste0(space_str, " ... +", length(issues$spaces) - 3, " more")
            }
            cat("    -", var_name, ": spaces in values -", space_str, "\n")
          }

          # Report underscores (can conflict with FAMD naming)
          if (!is.null(issues$underscores)) {
            underscore_examples <- head(issues$underscores, 3)
            underscore_str <- paste0("\"", underscore_examples, "\"", collapse = ", ")
            if (length(issues$underscores) > 3) {
              underscore_str <- paste0(underscore_str, " ... +", length(issues$underscores) - 3, " more")
            }
            cat("    -", var_name, ": underscores in values -", underscore_str, "\n")
          }

          # Report other special chars
          if (!is.null(issues$special)) {
            special_examples <- head(issues$special, 3)
            special_str <- paste0("\"", special_examples, "\"", collapse = ", ")
            if (length(issues$special) > 3) {
              special_str <- paste0(special_str, " ... +", length(issues$special) - 3, " more")
            }
            cat("    -", var_name, ": special characters -", special_str, "\n")
          }

          n_shown <- n_shown + 1
        }
        cat("\n")
        cat("  TIP: Use gsub() to replace spaces/special chars before FAMD.\n")
        cat("       Example: data$col <- gsub(\" \", \"_\", data$col)\n\n")
      }
    }
  }

  invisible(x)
}


#' Convert Columns to Factors Based on QC Report
#'
#' @description Helper function to convert potential categorical columns
#' (identified by HADES_qc_mixed_data) to factors. Can also be used to convert
#' any specified columns.
#'
#' @param data Data frame to modify.
#' @param qc_report Optional. A hades_qc object from HADES_qc_mixed_data().
#'   If provided, will use the potential_categorical list by default.
#' @param columns Character vector. Column names to convert. If NULL and
#'   qc_report is provided, uses all potential_categorical columns.
#'   If both qc_report and columns are provided, columns takes precedence.
#' @param ordered Logical or named list. If TRUE, creates ordered factors
#'   (for ordinal data). If a named list, specifies which columns should be
#'   ordered (e.g., list(severity = TRUE, group = FALSE)). Default = FALSE.
#' @param exclude Character vector. Column names to exclude from conversion
#'   even if flagged. Useful when you've reviewed the QC and determined some
#'   flagged columns are actually continuous.
#' @param verbose Logical. Print information about conversions. Default = TRUE.
#'
#' @return The data frame with specified columns converted to factors.
#'
#' @details
#' This function is intended to be used after reviewing the QC report from
#' HADES_qc_mixed_data(). The typical workflow is:
#' 1. Run HADES_qc_mixed_data() and review potential_categorical columns
#' 2. Decide which ones are truly categorical
#' 3. Use this function to convert them, excluding any that are actually continuous
#'
#' For ordinal data (e.g., severity scores 0-3), set ordered = TRUE to create
#' ordered factors. The levels will be sorted numerically for numeric columns.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert all flagged potential categorical columns
#' data_clean <- HADES_convert_to_factor(my_data, qc_report)
#'
#' # Convert specific columns only
#' data_clean <- HADES_convert_to_factor(my_data, columns = c("severity", "stage"))
#'
#' # Exclude some flagged columns (they're actually continuous)
#' data_clean <- HADES_convert_to_factor(my_data, qc_report, exclude = c("age_group"))
#'
#' # Create ordered factors for ordinal data
#' data_clean <- HADES_convert_to_factor(my_data, qc_report, ordered = TRUE)
#'
#' }
HADES_convert_to_factor <- function(data,
                                     qc_report = NULL,
                                     columns = NULL,
                                     ordered = FALSE,
                                     exclude = NULL,
                                     verbose = TRUE) {

  # Input validation

  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }

  # Determine which columns to convert
  if (is.null(columns)) {
    if (is.null(qc_report)) {
      stop("Must provide either qc_report or columns argument")
    }
    if (!inherits(qc_report, "hades_qc")) {
      stop("qc_report must be a hades_qc object from HADES_qc_mixed_data()")
    }
    if (is.null(qc_report$type_issues$potential_categorical) ||
        length(qc_report$type_issues$potential_categorical) == 0) {
      if (verbose) cat("[HADES] No potential categorical columns flagged in QC report.\n")
      return(data)
    }
    columns <- names(qc_report$type_issues$potential_categorical)
  }

  # Apply exclusions
  if (!is.null(exclude)) {
    columns <- setdiff(columns, exclude)
    if (verbose && length(exclude) > 0) {
      excluded_actual <- intersect(exclude, names(qc_report$type_issues$potential_categorical))
      if (length(excluded_actual) > 0) {
        cat("[HADES] Excluding:", paste(excluded_actual, collapse = ", "), "\n")
      }
    }
  }

  # Validate columns exist
  missing_cols <- setdiff(columns, names(data))
  if (length(missing_cols) > 0) {
    warning("Columns not found in data: ", paste(missing_cols, collapse = ", "))
    columns <- intersect(columns, names(data))
  }

  if (length(columns) == 0) {
    if (verbose) cat("[HADES] No columns to convert.\n")
    return(data)
  }

  if (verbose) {
    cat("[HADES] Converting", length(columns), "column(s) to factor:\n")
  }

  # Handle ordered parameter
  if (is.logical(ordered) && length(ordered) == 1) {
    # Single logical value - apply to all
    ordered_map <- setNames(rep(ordered, length(columns)), columns)
  } else if (is.list(ordered)) {
    # Named list - use specified values, default FALSE for unspecified
    ordered_map <- setNames(rep(FALSE, length(columns)), columns)
    for (col in names(ordered)) {
      if (col %in% columns) {
        ordered_map[col] <- ordered[[col]]
      }
    }
  } else {
    ordered_map <- setNames(rep(FALSE, length(columns)), columns)
  }

  # Convert each column
  for (col in columns) {
    col_data <- data[[col]]
    is_ordered <- ordered_map[col]

    # Get unique values for level ordering
    unique_vals <- unique(col_data[!is.na(col_data)])

    # Sort levels - numeric sort for numbers, alphabetic for others
    if (is.numeric(col_data)) {
      levels_sorted <- sort(unique_vals)
    } else {
      # Try numeric sort if values look numeric
      numeric_attempt <- suppressWarnings(as.numeric(as.character(unique_vals)))
      if (all(!is.na(numeric_attempt))) {
        levels_sorted <- unique_vals[order(numeric_attempt)]
      } else {
        levels_sorted <- sort(as.character(unique_vals))
      }
    }

    # Convert to factor
    data[[col]] <- factor(col_data, levels = levels_sorted, ordered = is_ordered)

    if (verbose) {
      ord_str <- if (is_ordered) " (ordered)" else ""
      cat("  -", col, ":", length(levels_sorted), "levels", ord_str, "\n")
    }
  }

  if (verbose) {
    cat("\n[HADES] Conversion complete.\n")
  }

  return(data)
}


#' Remove Variables with High Missing Data
#'
#' @description Simple helper to remove columns (variables) that exceed a
#' missing data threshold. This is a straightforward filter that doesn't
#' require a QC report.
#'
#' @param data Data frame to filter.
#' @param na_threshold Numeric (0-1). Columns with NA proportion exceeding
#'   this threshold will be removed. Default = 0.2 (20%).
#' @param protect Character vector. Column names that should never be removed,
#'   regardless of their NA proportion (e.g., sample ID column).
#' @param na_strings Character vector. Strings to treat as NA values (e.g.,
#'   "N/A", "NA", "None"). These will be counted as missing AND converted to
#'   actual NA in the returned data. Set to NULL to disable. Default includes
#'   common patterns.
#' @param verbose Logical. Print information about removed columns. Default = TRUE.
#'
#' @return The data frame with high-NA columns removed (and NA-like strings
#'   converted to NA if na_strings is not NULL).
#'
#' @details
#' This is a simple filtering function for quick data cleaning. For more
#' comprehensive QC assessment, use HADES_qc_mixed_data() first.
#'
#' The function calculates the proportion of NA values in each column and
#' removes those exceeding the threshold. Protected columns are never removed.
#'
#' When na_strings is provided, matching strings in character/factor columns
#' are treated as missing data for threshold calculation and converted to
#' actual NA in the output.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Remove columns with >20% missing
#' data_clean <- HADES_drop_high_na_variables(my_data)
#'
#' # Stricter threshold (>10% missing)
#' data_clean <- HADES_drop_high_na_variables(my_data, na_threshold = 0.1)
#'
#' # Protect the ID column from removal
#' data_clean <- HADES_drop_high_na_variables(my_data, protect = "sample_id")
#'
#' # Disable NA-string conversion
#' data_clean <- HADES_drop_high_na_variables(my_data, na_strings = NULL)
#'
#' }
HADES_drop_high_na_variables <- function(data,
                                          na_threshold = 0.2,
                                          protect = NULL,
                                          na_strings = c(
                                            "NA", "N/A", "na", "n/a",
                                            "None", "none", "NONE",
                                            "NULL", "null",
                                            "NaN", "nan",
                                            "", " ",
                                            ".", "-", "--",
                                            "missing", "Missing", "MISSING",
                                            "#N/A", "#NA", "#VALUE!"
                                          ),
                                          verbose = TRUE) {

  # Input validation
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }

  if (na_threshold < 0 || na_threshold > 1) {
    stop("na_threshold must be between 0 and 1")
  }

  if (ncol(data) == 0) {
    if (verbose) cat("[HADES] Data has no columns.\n")
    return(data)
  }

  # Convert NA-like strings to actual NA
  n_converted <- 0
  if (!is.null(na_strings) && length(na_strings) > 0) {
    for (col_name in names(data)) {
      col <- data[[col_name]]

      # Only process character or factor columns
      if (is.character(col) || is.factor(col)) {
        col_char <- as.character(col)
        col_trimmed <- trimws(col_char)

        # Find matches
        is_na_string <- !is.na(col_char) & col_trimmed %in% na_strings
        n_matches <- sum(is_na_string)

        if (n_matches > 0) {
          n_converted <- n_converted + n_matches

          # Convert to NA
          if (is.factor(col)) {
            # For factors, need to handle carefully
            col_char[is_na_string] <- NA
            data[[col_name]] <- factor(col_char, levels = setdiff(levels(col), na_strings))
          } else {
            col[is_na_string] <- NA
            data[[col_name]] <- col
          }
        }
      }
    }

    if (n_converted > 0 && verbose) {
      cat("[HADES] Converted", n_converted, "NA-like string(s) to NA\n\n")
    }
  }

  # Calculate NA proportion per column
  na_props <- colMeans(is.na(data))

  # Identify columns to remove
  high_na_cols <- names(na_props)[na_props > na_threshold]

  # Apply protection
  if (!is.null(protect)) {
    protected_high_na <- intersect(high_na_cols, protect)
    if (length(protected_high_na) > 0 && verbose) {
      cat("[HADES] Protected columns with high NA (not removed):\n")
      for (col in protected_high_na) {
        cat("  -", col, ":", round(na_props[col] * 100, 1), "% NA\n")
      }
      cat("\n")
    }
    high_na_cols <- setdiff(high_na_cols, protect)
  }

  if (length(high_na_cols) == 0) {
    if (verbose) {
      cat("[HADES] No columns exceed", na_threshold * 100, "% missing threshold.\n")
      cat("[HADES] All", ncol(data), "columns retained.\n")
    }
    return(data)
  }

  # Report what will be removed
  if (verbose) {
    cat("[HADES] Removing", length(high_na_cols), "column(s) with >",
        na_threshold * 100, "% missing:\n", sep = "")
    # Sort by NA proportion (highest first)
    high_na_cols_sorted <- high_na_cols[order(-na_props[high_na_cols])]
    n_show <- min(20, length(high_na_cols_sorted))
    for (col in high_na_cols_sorted[1:n_show]) {
      cat("  -", col, ":", round(na_props[col] * 100, 1), "% NA\n")
    }
    if (length(high_na_cols_sorted) > 20) {
      cat("  ... and", length(high_na_cols_sorted) - 20, "more\n")
    }
    cat("\n")
  }

  # Remove columns
  data_filtered <- data[, !names(data) %in% high_na_cols, drop = FALSE]

  if (verbose) {
    cat("[HADES] Result:", ncol(data_filtered), "of", ncol(data), "columns retained.\n")
  }

  return(data_filtered)
}
