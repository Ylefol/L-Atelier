#' Export Processed Mass Spec Data as Per-Sample Files for the Sample Sheet
#'
#' @description Writes a processed \code{massspec_data} object's final
#' values out as one plain \code{id}/\code{value} file per sample -- the
#' same shape \code{ELEUTHIA_load_count_file()} already reads -- so
#' proteomics can register in the same canonical sample sheet as
#' RNA-seq/ChIP-seq/ATAC-seq, despite having no \code{format} value of its
#' own in \code{ELEUTHIA_validate_sample_sheet()} for a wide matrix.
#'
#' Run this on a \code{massspec_data} object \strong{after} all cross-sample
#' processing (\code{HADES_filter_massspec()},
#' \code{HADES_filter_massspec_proteins()}, \code{POSEIDON_impute_massspec()},
#' batch correction, etc.) -- this function only exports; it does not
#' filter, impute, or correct anything itself. Pooled/control samples
#' (any \code{SampleType != "SAMPLE"}) are never exported, since they
#' aren't biological samples the rest of GAIA would analyse -- they're
#' meant to be consumed entirely at the QC stage, before this runs.
#'
#' @param massspec_data A \code{massspec_data} object, normally post
#'   \code{POSEIDON_impute_massspec()}.
#' @param output_dir Character. Directory to write one file per sample
#'   into. Created if it doesn't exist.
#' @param group_col Character. Column in \code{massspec_data$sample_meta}
#'   to use as the sample sheet's \code{group} value. Must already exist --
#'   add it to \code{sample_meta} before calling if it isn't there yet
#'   (mirrors \code{ARTEMIS_normalize_counts()}'s \code{group_col}
#'   convention: GAIA expects grouping to already live in the metadata,
#'   not be invented by the function doing the work).
#' @param omics Character. Value for the returned rows' \code{omics}
#'   column. Default \code{"Proteomics"}.
#' @param id_col Character or NULL. Column in
#'   \code{massspec_data$protein_meta} to use as each file's identifier
#'   column. \code{NULL} (default) uses \code{rownames(massspec_data$wide)}
#'   directly (the protein ID column given to the loader as \code{id_col}).
#'   If duplicated, values are
#'   made unique with \code{make.unique()} and a warning is issued.
#' @param value_col Character. Header for the value column in each written
#'   file. Default \code{"intensity"} -- a label only; reading back via
#'   \code{ELEUTHIA_load_count_file()} is positional, not name-matched.
#' @param file_suffix Character. Appended to each sample_id for the output
#'   filename: \code{<sample_id><file_suffix>}. Default
#'   \code{"_proteomics.txt"}.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#'
#' @return A data.frame of new sample sheet rows (\code{sample_id,
#'   file_loc, file_name, group, omics, format, bed_loc}; \code{format}
#'   is always \code{"counts"}, \code{bed_loc} always \code{NA}), ready to
#'   \code{rbind()} onto an existing sample sheet and re-validate with
#'   \code{ELEUTHIA_validate_sample_sheet()}. Does not write or modify any
#'   existing sample sheet file -- merging is the caller's decision.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ms <- HADES_filter_massspec(ms_raw)
#' ms <- HADES_filter_massspec_proteins(ms)
#' ms <- POSEIDON_impute_massspec(ms)
#' ms$sample_meta$group <- ifelse(ms$sample_meta$SampleID == "ACH-000551", "K562", "Other")
#'
#' new_rows <- ELEUTHIA_export_massspec_to_sheet(ms, "prepped/proteomics_counts",
#'                                                group_col = "group")
#'
#' sheet <- read.csv("sample_sheet.csv", stringsAsFactors = FALSE)
#' sheet <- rbind(sheet, new_rows)
#' write.csv(sheet, "sample_sheet.csv", row.names = FALSE)
#' }
ELEUTHIA_export_massspec_to_sheet <- function(massspec_data,
                                               output_dir,
                                               group_col,
                                               omics       = "Proteomics",
                                               id_col      = NULL,
                                               value_col   = "intensity",
                                               file_suffix = "_proteomics.txt",
                                               verbose     = TRUE) {

  if (!inherits(massspec_data, "massspec_data"))
    stop("massspec_data must be a massspec_data object.", call. = FALSE)

  smeta   <- massspec_data$sample_meta
  sid_col <- massspec_data$params$sample_col
  if (is.null(sid_col) || !sid_col %in% colnames(smeta))
    stop("massspec_data$params$sample_col not found in sample_meta.", call. = FALSE)

  if (!group_col %in% colnames(smeta))
    stop("group_col '", group_col, "' not found in massspec_data$sample_meta. ",
         "Add it before calling this function.", call. = FALSE)

  if ("SampleType" %in% colnames(smeta)) {
    keep <- smeta$SampleType == "SAMPLE"
    if (verbose && any(!keep))
      cat("[ELEUTHIA] Excluding", sum(!keep),
          "non-SAMPLE row(s) (pools/controls) from export\n")
    smeta <- smeta[keep, , drop = FALSE]
  }

  if (nrow(smeta) == 0)
    stop("No SAMPLE-type rows left to export.", call. = FALSE)

  wide  <- massspec_data$wide
  pmeta <- massspec_data$protein_meta

  ids <- if (is.null(id_col)) {
    rownames(wide)
  } else {
    if (!id_col %in% colnames(pmeta))
      stop("id_col '", id_col, "' not found in massspec_data$protein_meta.", call. = FALSE)
    id_vals <- pmeta[rownames(wide), id_col]
    if (anyDuplicated(id_vals)) {
      n_dup <- sum(duplicated(id_vals))
      warning(n_dup, " duplicate value(s) in id_col '", id_col, "' -- ",
              "made unique with make.unique(). Consider id_col = NULL to use ",
              "the loader's own (unique) row identifiers instead.", call. = FALSE)
      id_vals <- make.unique(as.character(id_vals), sep = "_dup")
    }
    id_vals
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (verbose)
    cat("[ELEUTHIA] Exporting", nrow(smeta), "sample(s) to", output_dir, "\n")

  sheet_rows <- vector("list", nrow(smeta))

  for (i in seq_len(nrow(smeta))) {
    sid <- as.character(smeta[[sid_col]][i])
    if (!sid %in% colnames(wide)) {
      warning("Sample '", sid, "' not found in massspec_data$wide -- skipped.", call. = FALSE)
      next
    }

    out_df <- data.frame(id = ids, value = wide[, sid])
    colnames(out_df) <- c("id", value_col)

    file_name <- paste0(sid, file_suffix)
    write.table(out_df, file.path(output_dir, file_name),
                sep = "\t", quote = FALSE, row.names = FALSE)

    sheet_rows[[i]] <- data.frame(
      sample_id = sid,
      file_loc  = output_dir,
      file_name = file_name,
      group     = as.character(smeta[[group_col]][i]),
      omics     = omics,
      format    = "counts",
      bed_loc   = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  result <- do.call(rbind, sheet_rows)
  if (is.null(result))
    stop("No samples were exported: none of the SAMPLE-type samples in ",
         "massspec_data$sample_meta were found in massspec_data$wide.", call. = FALSE)
  rownames(result) <- NULL

  if (verbose)
    cat("[ELEUTHIA] Wrote", nrow(result), "file(s); returning", nrow(result),
        "sample sheet row(s) (not written to disk -- rbind onto your sheet yourself)\n")

  result
}
