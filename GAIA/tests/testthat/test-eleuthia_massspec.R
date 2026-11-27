# ELEUTHIA mass spectrometry: the generic loader (a protein x sample table that has
# been prepared into an explicit shape) and the export of a processed object as
# per-sample files for the sample sheet. Data are simulated; values are chosen so the
# log2 transform and the exports have known answers.

# 6 proteins x (4 samples + 1 pool), with a summary row that has a blank protein ID.
# Sample columns follow a "{SampleID}_Plate N" naming, the way some pipelines name
# them; intensities are powers of two so log2() gives whole numbers.
make_pg_matrix <- function(dir, name = "report.pg_matrix.tsv", sep = "\t") {
  df <- data.frame(
    Protein.Group = c("P1", "P2", "P3", "P4", "P5", "P6", ""),
    Protein.Names = c("A_HUMAN", "B_HUMAN", "C_HUMAN", "D_HUMAN", "E_HUMAN", "F_HUMAN", ""),
    Genes         = c("A", "B", "C", "D", "E", "F", ""),
    First.Protein.Description = paste("protein", c("A", "B", "C", "D", "E", "F", "summary")),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  df[["1001_Plate 1"]]     <- c(1024, 2048, 4096, 0, 8, 1e6, 99)
  df[["1002_Plate 1"]]     <- c(1024, 2048, 4096, 16, NA, 1e6, 99)
  df[["SEP001-3_Plate 2"]] <- c(512, 1024, 2048, 32, 16, 1e6, 99)
  df[["SEP002-1_Plate 2"]] <- c(512, 1024, 2048, 64, 32, 1e6, 99)
  df[["pool1_Plate 1"]]    <- c(1000, 2000, 4000, 100, 10, 1e6, 99)
  path <- file.path(dir, name)
  utils::write.table(df, path, sep = sep, quote = FALSE, row.names = FALSE, na = "")
  path
}

make_pg_metadata <- function(dir, ids = c("1001", "1002", "SEP001-3")) {
  path <- file.path(dir, "meta.csv")
  utils::write.csv(data.frame(SampleID = ids,
                              Group = rep(c("case", "control"), length.out = length(ids)),
                              stringsAsFactors = FALSE),
                   path, row.names = FALSE)
  path
}

PLATE_RE <- "^(.+)_Plate [0-9]+$"

# the plate-named table through the generic loader (blank-ID summary row is dropped)
load_plate_ms <- function(dir, ..., scale = "linear", name = "report.pg_matrix.tsv") {
  quiet(ELEUTHIA_load_massspec(make_pg_matrix(dir, name), scale = scale,
                               sample_regex = PLATE_RE, pool_pattern = "^pool", ...))
}

# a small prepared matrix: ID column, two annotation columns, three samples
prepared_df <- function() {
  data.frame(protein = c("P1", "P2", "P3"), gene = c("A", "B", "C"), mw = c(10, 20, 30),
             S1 = c(1024, 2048, 0), S2 = c(2048, 4096, 8), S3 = c(512, NA, 16),
             stringsAsFactors = FALSE)
}

# ---- scale: required, explicit -----------------------------------------------

test_that("scale is required and must be log2 or linear", {
  expect_error(ELEUTHIA_load_massspec(prepared_df()), "'scale' is required")
  expect_error(ELEUTHIA_load_massspec(prepared_df(), scale = "ln"), "must be")
  expect_error(ELEUTHIA_load_massspec(prepared_df(), scale = c("log2", "linear")), "must be")
})

test_that("linear scale is log2-transformed, with zeros and negatives set to NA", {
  df <- prepared_df(); df$S3[1] <- -5
  ms <- quiet(ELEUTHIA_load_massspec(df, scale = "linear", annotation_cols = c("gene", "mw")))

  expect_equal(unname(ms$wide["P1", "S1"]), 10)          # log2(1024)
  expect_equal(unname(ms$wide["P2", "S2"]), 12)          # log2(4096)
  expect_equal(unname(ms$wide["P3", "S3"]), 4)           # log2(16)
  expect_true(is.na(ms$wide["P3", "S1"]))                # 0 -> NA, not -Inf
  expect_true(is.na(ms$wide["P1", "S3"]))                # negative -> NA
  expect_true(is.na(ms$wide["P2", "S3"]))                # missing stays NA
  expect_false(any(is.infinite(ms$wide)))
  expect_equal(ms$params$n_nonpositive_values_set_na, 2L)
  expect_equal(ms$params$scale, "linear")
})

test_that("log2 scale keeps the values as given, including negatives", {
  df <- data.frame(protein = c("P1", "P2"), S1 = c(20.5, -1.2), S2 = c(NA, 18))
  ms <- quiet(ELEUTHIA_load_massspec(df, scale = "log2"))

  expect_equal(unname(ms$wide["P1", "S1"]), 20.5)
  expect_equal(unname(ms$wide["P2", "S1"]), -1.2)
  expect_true(is.na(ms$wide["P1", "S2"]))
})

test_that("a scale that looks wrong for the values raises a warning", {
  lin <- data.frame(protein = c("P1", "P2"), S1 = c(1e6, 2e6), S2 = c(3e6, 4e6))
  expect_warning(quiet(ELEUTHIA_load_massspec(lin, scale = "log2")), "linear scale")

  logged <- data.frame(protein = c("P1", "P2"), S1 = c(20, 22), S2 = c(21, 25))
  expect_warning(quiet(ELEUTHIA_load_massspec(logged, scale = "linear")), "already be log")
})

# ---- shape and column roles ---------------------------------------------------

test_that("a prepared matrix loads: ID column, annotation columns, samples", {
  ms <- quiet(ELEUTHIA_load_massspec(prepared_df(), scale = "linear",
                                     annotation_cols = c("gene", "mw")))

  expect_s3_class(ms, "massspec_data")
  expect_named(ms, c("wide", "sample_meta", "protein_meta", "params"))
  expect_equal(dim(ms$wide), c(3L, 3L))
  expect_equal(rownames(ms$wide), c("P1", "P2", "P3"))
  expect_equal(colnames(ms$wide), c("S1", "S2", "S3"))
  expect_equal(colnames(ms$protein_meta), c("protein", "gene", "mw"))
  expect_equal(rownames(ms$protein_meta), rownames(ms$wide))
  expect_equal(ms$protein_meta$gene, c("A", "B", "C"))
})

test_that("without sample_regex every column except the ID is a sample", {
  df <- data.frame(protein = c("P1", "P2"), S1 = c(1024, 2048), S2 = c(4096, 8192))
  ms <- quiet(ELEUTHIA_load_massspec(df, scale = "linear"))

  expect_equal(colnames(ms$wide), c("S1", "S2"))
  expect_equal(colnames(ms$protein_meta), "protein")
})

test_that("a non-numeric column that is not declared as annotation is an error naming it", {
  expect_error(ELEUTHIA_load_massspec(prepared_df(), scale = "linear"), "not numeric")
  expect_error(ELEUTHIA_load_massspec(prepared_df(), scale = "linear"), "gene")
})

test_that("sample_regex selects the sample columns and a capture group names them", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  ms <- load_plate_ms(d)

  expect_equal(dim(ms$wide), c(6L, 5L))                    # the blank-ID summary row is dropped
  expect_equal(rownames(ms$wide), paste0("P", 1:6))
  expect_equal(colnames(ms$wide), c("1001", "1002", "SEP001-3", "SEP002-1", "pool1"))
  # everything that is not a sample column is annotation
  expect_equal(colnames(ms$protein_meta),
               c("Protein.Group", "Protein.Names", "Genes", "First.Protein.Description"))
  expect_equal(ms$params$n_blank_id_rows_dropped, 1L)
})

test_that("the original column names are kept in source_column", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  sm <- load_plate_ms(d)$sample_meta

  expect_equal(sm$source_column, c("1001_Plate 1", "1002_Plate 1", "SEP001-3_Plate 2",
                                   "SEP002-1_Plate 2", "pool1_Plate 1"))
  # information in the name can be extracted afterwards, e.g. the plate
  plate <- sub(".*_(Plate [0-9]+)$", "\\1", sm$source_column)
  expect_equal(plate, c("Plate 1", "Plate 1", "Plate 2", "Plate 2", "Plate 1"))
  expect_equal(rownames(sm), sm$SampleID)
})

test_that("a sample_regex without a capture group uses the whole column name", {
  df <- data.frame(protein = c("P1", "P2"), a_int = c(1024, 2048), b_int = c(4096, 8192), note = c(1, 2))
  ms <- quiet(ELEUTHIA_load_massspec(df, scale = "linear", sample_regex = "_int$"))

  expect_equal(colnames(ms$wide), c("a_int", "b_int"))
  expect_equal(colnames(ms$protein_meta), c("protein", "note"))
})

test_that("MaxQuant-style columns: only the chosen intensity family becomes samples", {
  pg <- data.frame(`Protein IDs` = c("P1", "P2", "P3"), `Gene names` = c("A", "B", "C"),
                   `Potential contaminant` = c("", "+", ""),
                   `LFQ intensity 1_X` = c(1024, 2048, 0), `LFQ intensity 2_X` = c(2048, 4096, 16),
                   `Intensity 1_X` = c(1, 2, 3), `Intensity 2_X` = c(4, 5, 6),
                   check.names = FALSE, stringsAsFactors = FALSE)
  pg <- pg[!(pg$`Potential contaminant` %in% "+"), ]      # the user drops flagged rows first
  ms <- quiet(ELEUTHIA_load_massspec(pg, scale = "linear", id_col = "Protein IDs",
                                     sample_regex = "^LFQ intensity (.+)$"))

  expect_equal(colnames(ms$wide), c("1_X", "2_X"))
  expect_equal(rownames(ms$wide), c("P1", "P3"))
  expect_equal(unname(ms$wide["P1", "1_X"]), 10)
  expect_true(all(c("Gene names", "Intensity 1_X") %in% colnames(ms$protein_meta)))
})

test_that("sample_regex that matches nothing, and the regex + annotation_cols mix, are errors", {
  expect_error(ELEUTHIA_load_massspec(prepared_df(), scale = "linear", sample_regex = "^zzz"),
               "matched none")
  expect_error(ELEUTHIA_load_massspec(prepared_df(), scale = "linear",
                                      sample_regex = "^S", annotation_cols = "gene"),
               "not both")
})

test_that("id_col can be a name or a position, and a bad one is an error", {
  df <- data.frame(S1 = c(1024, 2048), acc = c("P1", "P2"), S2 = c(4096, 8192))
  by_name <- quiet(ELEUTHIA_load_massspec(df, scale = "linear", id_col = "acc"))
  by_pos  <- quiet(ELEUTHIA_load_massspec(df, scale = "linear", id_col = 2))

  expect_equal(rownames(by_name$wide), c("P1", "P2"))
  expect_equal(by_pos$wide, by_name$wide)
  expect_equal(colnames(by_name$wide), c("S1", "S2"))
  expect_error(ELEUTHIA_load_massspec(df, scale = "linear", id_col = "nope"), "not found")
  expect_error(ELEUTHIA_load_massspec(df, scale = "linear", id_col = 9), "outside")
})

# ---- strictness: duplicates and bad values -----------------------------------

test_that("duplicated protein IDs are an error that names them", {
  df <- data.frame(protein = c("P1", "P1", "P2"), S1 = c(8, 16, 32))
  expect_error(ELEUTHIA_load_massspec(df, scale = "linear"), "Duplicated protein identifier")
  expect_error(ELEUTHIA_load_massspec(df, scale = "linear"), "P1")
})

test_that("duplicated sample IDs are an error that names them", {
  df <- data.frame(protein = c("P1", "P2"), `a_1` = c(1024, 2048), `a_2` = c(4096, 8192), check.names = FALSE)
  expect_error(ELEUTHIA_load_massspec(df, scale = "linear", sample_regex = "^(a)_[0-9]$"),
               "Duplicated sample ID")
})

test_that("recognised missing-value strings are read as NA and other text is an error", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  p <- file.path(d, "t.tsv")
  writeLines(c("protein\tS1\tS2", "P1\t1024\tFiltered", "P2\tNaN\t#N/A", "P3\t2048\t4096"), p)
  ms <- quiet(ELEUTHIA_load_massspec(p, scale = "linear"))
  expect_equal(unname(ms$wide["P1", "S1"]), 10)
  expect_equal(sum(is.na(ms$wide)), 3L)

  writeLines(c("protein\tS1\tS2", "P1\t1024\tlow", "P2\t2048\t4096"), p)
  expect_error(ELEUTHIA_load_massspec(p, scale = "linear"), "not numeric")
})

# ---- input forms --------------------------------------------------------------

test_that("csv, tsv and txt files and a data.frame give the same object", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  df <- prepared_df()
  utils::write.csv(df, file.path(d, "a.csv"), row.names = FALSE)
  utils::write.table(df, file.path(d, "a.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(df, file.path(d, "a.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  load <- function(x) quiet(ELEUTHIA_load_massspec(x, scale = "linear", annotation_cols = c("gene", "mw")))
  ref <- load(df)

  for (f in c("a.csv", "a.tsv", "a.txt")) {
    got <- load(file.path(d, f))
    expect_equal(got$wide, ref$wide, info = f)
    expect_equal(got$protein_meta, ref$protein_meta, info = f)
  }
})

test_that("bad inputs are rejected with a clear message", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  bad <- file.path(d, "x.parquet"); writeLines("x", bad)

  expect_error(ELEUTHIA_load_massspec(file.path(d, "nope.tsv"), scale = "linear"), "not found")
  expect_error(ELEUTHIA_load_massspec(bad, scale = "linear"), "Unsupported format")
  expect_error(ELEUTHIA_load_massspec(42, scale = "linear"), "path or a data.frame")
  expect_error(ELEUTHIA_load_massspec(data.frame(protein = c("", NA), S1 = c(1, 2)), scale = "linear"),
               "No rows left")
})

# ---- sample types and metadata -----------------------------------------------

test_that("without metadata every sample is SAMPLE, and pool_pattern marks the pools", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  expect_equal(load_plate_ms(d)$sample_meta$SampleType,
               c("SAMPLE", "SAMPLE", "SAMPLE", "SAMPLE", "POOL"))

  no_pool <- quiet(ELEUTHIA_load_massspec(make_pg_matrix(d), scale = "linear", sample_regex = PLATE_RE))
  expect_true(all(no_pool$sample_meta$SampleType == "SAMPLE"))
})

test_that("metadata is joined; samples missing from it are OTHER; pools stay POOL", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  ms <- load_plate_ms(d, metadata = make_pg_metadata(d))
  sm <- ms$sample_meta

  expect_true("Group" %in% colnames(sm))
  expect_equal(sm["1001", "SampleType"], "SAMPLE")
  expect_equal(sm["SEP001-3", "SampleType"], "SAMPLE")
  expect_equal(sm["SEP002-1", "SampleType"], "OTHER")      # not in the metadata
  expect_equal(sm["pool1", "SampleType"], "POOL")          # pools never need metadata
  expect_equal(sm["1001", "Group"], "case")
  expect_true(is.na(sm["SEP002-1", "Group"]))
  expect_equal(ncol(ms$wide), 5L)                          # the matrix is unaffected
})

test_that("a SampleType column in the metadata is used for the samples it lists", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  meta <- data.frame(SampleID = c("1001", "1002", "SEP001-3", "SEP002-1"),
                     SampleType = c("SAMPLE", "SAMPLE", "REFERENCE", "SAMPLE"))
  sm <- load_plate_ms(d, metadata = meta)$sample_meta

  expect_equal(sm["SEP001-3", "SampleType"], "REFERENCE")
  expect_equal(sm["pool1", "SampleType"], "POOL")
  expect_equal(sum(colnames(sm) == "SampleType"), 1L)      # joined once, not duplicated
})

test_that("a missing SampleType for a listed sample, or duplicate metadata IDs, are errors", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  meta <- data.frame(SampleID = c("1001", "1002"), SampleType = c("SAMPLE", NA))
  expect_error(load_plate_ms(d, metadata = meta), "Missing 'SampleType'")

  dup <- data.frame(SampleID = c("1001", "1001"), Group = c("a", "b"))
  expect_error(load_plate_ms(d, metadata = dup), "Duplicated sample ID")
})

test_that("a custom sample_col names the sample ID column, in the metadata and the result", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  meta <- data.frame(Patient = c("1001", "1002"), Arm = c("x", "y"))
  ms <- load_plate_ms(d, metadata = meta, sample_col = "Patient")

  expect_true("Patient" %in% colnames(ms$sample_meta))
  expect_equal(ms$sample_meta["1002", "Arm"], "y")
  expect_equal(ms$params$sample_col, "Patient")
  expect_error(load_plate_ms(d, metadata = meta, sample_col = "SampleID"), "not found")
})

test_that("metadata can be a path or a data.frame with the same result", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  p <- make_pg_metadata(d)
  from_path <- load_plate_ms(d, metadata = p)
  from_df   <- load_plate_ms(d, metadata = utils::read.csv(p, stringsAsFactors = FALSE))

  expect_equal(from_path$sample_meta, from_df$sample_meta)
})

test_that("params record how the object was built", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  ms <- load_plate_ms(d)

  expect_equal(ms$params$id_col, "Protein.Group")
  expect_equal(ms$params$sample_regex, PLATE_RE)
  expect_equal(ms$params$pool_pattern, "^pool")
  expect_false("source" %in% names(ms$params))
  expect_equal(ms$params$ms_data, file.path(d, "report.pg_matrix.tsv"))
  expect_equal(ms$params$metadata, NULL)
})

test_that("print.massspec_data summarises the object", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  ms <- load_plate_ms(d)
  out <- utils::capture.output(res <- print(ms))

  expect_identical(res, ms)       # returns its input invisibly
  expect_true(any(grepl("Proteins : 6", out)))
  expect_true(any(grepl("Samples  : 5", out)))
  expect_true(any(grepl("SAMPLE=4", out)))
  expect_false(any(grepl("DIA", out)))
})

test_that("the loader has no source argument and no tool-specific arguments", {
  expect_false("source" %in% names(formals(ELEUTHIA_load_massspec)))
  expect_false("n_annotation_cols" %in% names(formals(ELEUTHIA_load_massspec)))
  expect_false(any(c("ms_file", "metadata_file") %in% names(formals(ELEUTHIA_load_massspec))))
})

# ---- export to the sample sheet -----------------------------------------------

# A processed object: group lives in sample_meta before export.
make_ms_for_export <- function() {
  d <- tempfile("ms_src_"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  ms <- load_plate_ms(d, metadata = make_pg_metadata(d, ids = c("1001", "1002", "SEP001-3", "SEP002-1")))
  ms$sample_meta$group <- c("case", "case", "control", "control", NA)
  ms
}

export_ms <- function(ms, out, ...) {
  quiet(ELEUTHIA_export_massspec_to_sheet(ms, out, group_col = "group", ...))
}

test_that("export writes one file per SAMPLE and never exports pools", {
  ms <- make_ms_for_export()
  out <- file.path(tempfile("ms_out_"), "nested"); on.exit(unlink(dirname(out), recursive = TRUE), add = TRUE)
  rows <- export_ms(ms, out)

  expect_true(dir.exists(out))                    # created, including the parent
  expect_equal(sort(list.files(out)),
               sort(paste0(c("1001", "1002", "SEP001-3", "SEP002-1"), "_proteomics.txt")))
  expect_false(any(grepl("pool", list.files(out))))
  expect_equal(nrow(rows), 4L)
})

test_that("the returned rows are a valid sample sheet block", {
  ms <- make_ms_for_export()
  out <- tempfile("ms_out_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  rows <- export_ms(ms, out)

  expect_named(rows, c("sample_id", "file_loc", "file_name", "group", "omics", "format", "bed_loc"))
  expect_equal(rows$sample_id, c("1001", "1002", "SEP001-3", "SEP002-1"))
  expect_equal(rows$group, c("case", "case", "control", "control"))
  expect_true(all(rows$omics == "Proteomics"))
  expect_true(all(rows$format == "counts"))
  expect_true(all(is.na(rows$bed_loc)))
  expect_equal(rownames(rows), as.character(1:4))
  expect_true(all(file.exists(file.path(rows$file_loc, rows$file_name))))

  v <- quiet(ELEUTHIA_validate_sample_sheet(rows, check_beds = FALSE, verbose = FALSE))
  expect_true(v$valid)
})

test_that("an exported file reads back through load_count_file with the same values", {
  ms <- make_ms_for_export()
  out <- tempfile("ms_out_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  rows <- export_ms(ms, out)

  back <- ELEUTHIA_load_count_file(file.path(out, rows$file_name[rows$sample_id == "1001"]))
  expect_equal(names(back), rownames(ms$wide))
  expect_equal(unname(back), unname(ms$wide[, "1001"]))
  expect_true(is.na(back["P4"]))                  # NA survives the round trip

  hdr <- readLines(file.path(out, "1001_proteomics.txt"), n = 1)
  expect_equal(hdr, "id\tintensity")
})

test_that("the exported files load as one matrix via the sample sheet", {
  ms <- make_ms_for_export()
  out <- tempfile("ms_out_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  rows <- export_ms(ms, out)

  loaded <- quiet(ELEUTHIA_load_rnaseq_from_sheet(rows, omics = "Proteomics", verbose = FALSE))
  expect_equal(dim(loaded$counts), c(6L, 4L))
  expect_equal(colnames(loaded$counts), rows$sample_id)
  expect_equal(unname(loaded$counts[, "SEP001-3"]), unname(ms$wide[, "SEP001-3"]))
})

test_that("omics, value_col and file_suffix control the output naming", {
  ms <- make_ms_for_export()
  out <- tempfile("ms_out_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  rows <- export_ms(ms, out, omics = "Phospho", value_col = "log2_int", file_suffix = ".tsv")

  expect_true(all(rows$omics == "Phospho"))
  expect_equal(rows$file_name[1], "1001.tsv")
  expect_equal(readLines(file.path(out, "1001.tsv"), n = 1), "id\tlog2_int")
})

test_that("id_col swaps the identifier column, and duplicates are made unique with a warning", {
  ms <- make_ms_for_export()
  out <- tempfile("ms_out_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)

  rows <- export_ms(ms, out, id_col = "Genes")
  back <- ELEUTHIA_load_count_file(file.path(out, rows$file_name[1]))
  expect_equal(names(back), c("A", "B", "C", "D", "E", "F"))

  ms$protein_meta$Genes[2] <- "A"       # now A appears twice
  out2 <- tempfile("ms_out2_"); on.exit(unlink(out2, recursive = TRUE), add = TRUE)
  expect_warning(export_ms(ms, out2, id_col = "Genes"), "duplicate")
  back2 <- ELEUTHIA_load_count_file(file.path(out2, "1001_proteomics.txt"))
  expect_false(anyDuplicated(names(back2)) > 0)
  expect_equal(names(back2)[2], "A_dup1")
})

test_that("the export does not modify the object or write a sample sheet", {
  ms <- make_ms_for_export()
  before <- ms
  out <- tempfile("ms_out_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  export_ms(ms, out)

  expect_identical(ms, before)
  expect_false(any(grepl("\\.csv$", list.files(out))))
})

test_that("a sample missing from the matrix is skipped with a warning", {
  ms <- make_ms_for_export()
  ms$wide <- ms$wide[, colnames(ms$wide) != "1002", drop = FALSE]
  out <- tempfile("ms_out_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)

  expect_warning(rows <- export_ms(ms, out), "1002")
  expect_equal(rows$sample_id, c("1001", "SEP001-3", "SEP002-1"))
  expect_false(file.exists(file.path(out, "1002_proteomics.txt")))
})

test_that("export validates its inputs", {
  ms <- make_ms_for_export()
  out <- tempfile("ms_out_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)

  expect_error(ELEUTHIA_export_massspec_to_sheet(list(), out, group_col = "group"),
               "massspec_data object")
  expect_error(ELEUTHIA_export_massspec_to_sheet(ms, out, group_col = "nope", verbose = FALSE),
               "group_col")
  expect_error(export_ms(ms, out, id_col = "nope"), "id_col")

  no_samples <- ms
  no_samples$sample_meta$SampleType <- "POOL"
  expect_error(export_ms(no_samples, out), "No SAMPLE-type rows")

  bad_col <- ms
  bad_col$params$sample_col <- "Missing"
  expect_error(export_ms(bad_col, out), "sample_col")
})

test_that("export stops with a clear error when none of the samples can be exported", {
  ms <- make_ms_for_export()
  ms$wide <- ms$wide[, c("pool1"), drop = FALSE]      # no SAMPLE column left in the matrix
  out <- tempfile("ms_out_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)

  expect_error(suppressWarnings(export_ms(ms, out)), "No samples were exported")
})
