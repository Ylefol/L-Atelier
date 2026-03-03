# Internal helper functions — not exported

# Extract a single sample row from the sample sheet.
# Stops with an informative message if the sample ID is not found.
.get_sample_row <- function(sample_sheet, sample_id) {
  row <- sample_sheet[sample_sheet$sample_id == sample_id, ]
  if (nrow(row) == 0) stop("Sample ID not found in sample sheet: '", sample_id, "'")
  if (nrow(row) > 1)  stop("Duplicate sample ID in sample sheet: '", sample_id, "'")
  row
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
