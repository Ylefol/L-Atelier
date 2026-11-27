# Simulated Olink NPX data with planted QC problems, written as a CSV so the real
# loader (ELEUTHIA_load_olink) is exercised without needing arrow/parquet.
#
# 30 proteins (OID001..OID030) x 25 samples:
#   S01-S22  biological samples (SampleType SAMPLE); S21, S22 have SampleQC FAIL
#   N1, N2   negative controls
#   X1       a sample that is not in the metadata file
# Planted problems (known answers):
#   OID001, OID002  NA in half of the biological samples      -> NA fraction 0.5
#   OID003          AssayQC WARN in 12 of 24 samples          -> warn fraction 0.5
#   OID004          AssayQC WARN in 1 sample                  -> kept at 0.2
#   S05             shifted +8 NPX on every protein           -> the outlier
# Plus 25 instrument-control rows (AssayType "ext_ctrl") and 3 rows with a blank
# SampleID, which the loader drops.
make_npx_long <- function(seed = 11) {
  set.seed(seed)
  samples  <- c(sprintf("S%02d", 1:22), "N1", "N2", "X1")
  proteins <- sprintf("OID%03d", 1:30)

  g <- expand.grid(SampleID = samples, OlinkID = proteins, stringsAsFactors = FALSE)
  g$Assay    <- paste0("GENE", sub("OID", "", g$OlinkID))
  g$UniProt  <- paste0("P", sub("OID", "", g$OlinkID))
  g$Panel    <- "Explore"
  g$NPX      <- stats::rnorm(nrow(g), mean = 5, sd = 1)
  g$AssayQC  <- "PASS"
  g$SampleQC <- ifelse(g$SampleID %in% c("S21", "S22"), "FAIL", "PASS")
  g$AssayType  <- "assay"
  g$SampleType <- ifelse(g$SampleID %in% c("N1", "N2"), "NEGATIVE_CONTROL", "SAMPLE")
  g$PlateID  <- ifelse(match(g$SampleID, samples) %% 2 == 0, "Plate1", "Plate2")
  g$WellID   <- paste0("A", match(g$SampleID, samples))

  g$NPX[g$SampleID == "S05"] <- g$NPX[g$SampleID == "S05"] + 8
  g$NPX[g$OlinkID %in% c("OID001", "OID002") & g$SampleID %in% sprintf("S%02d", 1:10)] <- NA
  g$AssayQC[g$OlinkID == "OID003" & g$SampleID %in% sprintf("S%02d", 1:12)] <- "WARN"
  g$AssayQC[g$OlinkID == "OID004" & g$SampleID == "S01"] <- "WARN"

  ctrl <- g[g$OlinkID == "OID001", ]
  ctrl$OlinkID <- "CTRL1"; ctrl$AssayType <- "ext_ctrl"
  blank <- g[g$OlinkID == "OID001", ][1:3, ]
  blank$SampleID <- c("", NA, "  ")
  rbind(g, ctrl, blank)
}

# Writes the NPX CSV and a metadata CSV (all samples except X1) into `dir`.
write_olink_files <- function(dir, seed = 11) {
  npx  <- file.path(dir, "npx.csv")
  meta <- file.path(dir, "meta.csv")
  utils::write.csv(make_npx_long(seed), npx, row.names = FALSE)
  ids <- c(sprintf("S%02d", 1:22), "N1", "N2")
  utils::write.csv(data.frame(SampleID = ids,
                              Group = rep(c("case", "control"), length.out = length(ids)),
                              stringsAsFactors = FALSE),
                   meta, row.names = FALSE)
  list(npx = npx, meta = meta)
}

# A freshly loaded object (restrict_to_metadata = TRUE drops X1).
load_sim_olink <- function(seed = 11) {
  d <- tempfile("olink_"); dir.create(d)
  on.exit(unlink(d, recursive = TRUE))
  f <- write_olink_files(d, seed)
  quiet(ELEUTHIA_load_olink(f$npx, metadata_file = f$meta,
                            restrict_to_metadata = TRUE, verbose = TRUE))
}

# The full chain used by most tests, with the outlier detected and removed.
run_olink_chain <- function(ol_raw) {
  ol <- quiet(HADES_filter_olink(ol_raw))
  ol <- quiet(HADES_filter_olink_proteins(ol, max_na_fraction = 0.2))
  ol <- quiet(HADES_filter_olink_assay_warn(ol, max_warn_fraction = 0.2))
  ol <- quiet(HADES_detect_outliers_olink(ol, n_pcs = 3, threshold = 3))
  quiet(HADES_detect_outliers_olink(ol, n_pcs = 3, threshold = 3, filter = TRUE))
}
