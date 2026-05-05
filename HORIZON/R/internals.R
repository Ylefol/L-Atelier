# Internal helper functions — not exported

# Extract a single sample row from the sample sheet.
# Stops with an informative message if the sample ID is not found.
.get_sample_row <- function(sample_sheet, sample_id) {
  row <- sample_sheet[sample_sheet$sample_id == sample_id, ]
  if (nrow(row) == 0) stop("Sample ID not found in sample sheet: '", sample_id, "'")
  if (nrow(row) > 1)  stop("Duplicate sample ID in sample sheet: '", sample_id, "'")
  row
}

# Validate chemistry and kit columns when present in a sample sheet.
# Called by both HORIZON_validate_sample_sheet() (opportunistic) and
# HORIZON_validate_parse_sheet() (required columns already enforced upstream).
# Stops with an error on any invalid value.
.validate_chemistry_kit <- function(ss) {
  valid_chemistry <- c("v2", "v3")
  bad_chem <- setdiff(unique(ss$chemistry[nzchar(ss$chemistry)]), valid_chemistry)
  if (length(bad_chem) > 0)
    stop("Invalid chemistry value(s): ", paste(bad_chem, collapse = ", "),
         ". Must be one of: ", paste(valid_chemistry, collapse = ", "),
         call. = FALSE)

  valid_kits <- c("WT", "WT_mini", "WT_mega", "WT_mega_384", "WT_penta", "WT_penta_384")
  bad_kit <- setdiff(unique(ss$kit[nzchar(ss$kit)]), valid_kits)
  if (length(bad_kit) > 0)
    stop("Invalid kit value(s): ", paste(bad_kit, collapse = ", "),
         ". Must be one of: ", paste(valid_kits, collapse = ", "),
         call. = FALSE)
}

# Run an external command with conda env PATH injection.
# When a HORIZON log is active, subprocess stdout+stderr are tee'd to the log
# file so they appear in both the console and the log simultaneously.
# Without an active log, falls back to plain system2().
.horizon_run_with_log <- function(bin, args, conda_env) {
  conda_path <- paste(file.path(conda_env, "bin"), Sys.getenv("PATH"), sep = ":")
  log_file   <- HORIZON_get_log_file()

  if (!is.null(log_file)) {
    cmd <- paste(
      "env", paste0("PATH=", shQuote(conda_path)),
      shQuote(bin), paste(vapply(args, shQuote, character(1L)), collapse = " "),
      "2>&1 | tee -a", shQuote(log_file)
    )
    system(cmd)
  } else {
    system2(bin, args = args, env = paste0("PATH=", conda_path))
  }
}

# Map human-readable strandedness labels to featureCounts integer codes.
# Returns 0 (unstranded), 1 (forward), or 2 (reverse).
.strand_to_int <- function(strandedness) {
  map <- c("unstranded" = 0L, "forward" = 1L, "reverse" = 2L)
  val <- map[strandedness]
  if (is.na(val)) {
    stop("Unknown strandedness value: '", strandedness,
         "'. Must be one of: unstranded, forward, reverse.")
  }
  val
}
