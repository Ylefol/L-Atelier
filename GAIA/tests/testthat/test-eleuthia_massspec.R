# ELEUTHIA mass spectrometry: loading a DIA-NN report.pg_matrix (with the
# `source` dispatch) and exporting a processed massspec_data object as
# per-sample files that register in the sample sheet. Data are simulated; the
# values are chosen so the log2 transform and the exports have known answers.

# 6 proteins x (4 samples + 1 pool). Column names follow the
# "{SampleID}_Plate N" convention the loader documents.
# Intensities are powers of two so log2() gives whole numbers.
make_pg_matrix <- function(dir, name = "report.pg_matrix.tsv", sep = "\t") {
  df <- data.frame(
    Protein.Group = c("P1", "P2", "P3", "P4", "P5", "P6", ""),
    Protein.Names = c("A_HUMAN", "B_HUMAN", "C_HUMAN", "D_HUMAN", "E_HUMAN", "F_HUMAN", ""),
    Genes         = c("A", "B", "C", "D", "E", "F", ""),
    First.Protein.Description = paste("protein", c("A", "B", "C", "D", "E", "F", "summary")),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  df[["1001_Plate 1"]]  <- c(1024, 2048, 4096, 0, 8, 1e6, 99)
  df[["1002_Plate 1"]]  <- c(1024, 2048, 4096, 16, NA, 1e6, 99)
  df[["SEP001-3_Plate 2"]] <- c(512, 1024, 2048, 32, 16, 1e6, 99)
  df[["SEP002-1_Plate 2"]] <- c(512, 1024, 2048, 64, 32, 1e6, 99)
  df[["pool1_Plate 1"]] <- c(1000, 2000, 4000, 100, 10, 1e6, 99)
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

load_ms <- function(...) quiet(ELEUTHIA_load_massspec(..., verbose = TRUE))

# ---- loading ------------------------------------------------------------------

test_that("load_massspec returns a massspec_data object with the right shape", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  ms <- load_ms(make_pg_matrix(d))

  expect_s3_class(ms, "massspec_data")
  expect_named(ms, c("wide", "sample_meta", "protein_meta", "params"))
  # 6 real proteins; the empty-Protein.Group summary row is dropped
  expect_equal(dim(ms$wide), c(6L, 5L))
  expect_equal(rownames(ms$wide), paste0("P", 1:6))
  expect_equal(colnames(ms$wide), c("1001", "1002", "SEP001-3", "SEP002-1", "pool1"))
  expect_equal(rownames(ms$protein_meta), rownames(ms$wide))
  expect_true(all(c("Protein.Group", "Genes") %in% colnames(ms$protein_meta)))
  expect_equal(ms$params$source, "dia-nn")
})

test_that("intensities are log2-transformed and zeros / missing values become NA", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  ms <- load_ms(make_pg_matrix(d))

  expect_equal(unname(ms$wide["P1", "1001"]), 10)         # log2(1024)
  expect_equal(unname(ms$wide["P3", "SEP001-3"]), 11)     # log2(2048)
  expect_equal(unname(ms$wide["P5", "1001"]), 3)          # log2(8)
  expect_true(is.na(ms$wide["P4", "1001"]))               # 0 -> NA, not -Inf
  expect_true(is.na(ms$wide["P5", "1002"]))               # missing stays NA
  expect_false(any(is.infinite(ms$wide)))
})

test_that("sample and plate IDs are parsed from the column names", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  ms <- load_ms(make_pg_matrix(d))
  sm <- ms$sample_meta

  expect_equal(sm$SampleID, c("1001", "1002", "SEP001-3", "SEP002-1", "pool1"))
  expect_equal(sm$PlateID, c("Plate 1", "Plate 1", "Plate 2", "Plate 2", "Plate 1"))
  expect_equal(rownames(sm), sm$SampleID)
})

test_that("without a metadata file every non-pool column is a SAMPLE and pools are POOL", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  sm <- load_ms(make_pg_matrix(d))$sample_meta

  expect_equal(sm$SampleType, c("SAMPLE", "SAMPLE", "SAMPLE", "SAMPLE", "POOL"))
})

test_that("a metadata file is joined, and unmatched non-pool columns become OTHER", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  ms <- load_ms(make_pg_matrix(d), metadata_file = make_pg_metadata(d))
  sm <- ms$sample_meta

  expect_true("Group" %in% colnames(sm))
  expect_equal(sm["1001", "SampleType"], "SAMPLE")
  expect_equal(sm["SEP001-3", "SampleType"], "SAMPLE")
  expect_equal(sm["SEP002-1", "SampleType"], "OTHER")     # not in the metadata
  expect_equal(sm["pool1", "SampleType"], "POOL")         # pools never need metadata
  expect_equal(sm["1001", "Group"], "case")
  expect_true(is.na(sm["SEP002-1", "Group"]))
  # the intensity matrix is unaffected by the join
  expect_equal(ncol(ms$wide), 5L)
})

test_that("a custom sample_col names the sample ID column", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  meta <- file.path(d, "meta2.csv")
  utils::write.csv(data.frame(Patient = c("1001", "1002"), Arm = c("x", "y")),
                   meta, row.names = FALSE)
  ms <- load_ms(make_pg_matrix(d), metadata_file = meta, sample_col = "Patient")

  expect_true("Patient" %in% colnames(ms$sample_meta))
  expect_equal(ms$sample_meta["1002", "Arm"], "y")
  expect_equal(ms$params$sample_col, "Patient")
})

test_that("tab, csv and txt files all load to the same matrix", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  tsv <- load_ms(make_pg_matrix(d, "a.tsv"))
  txt <- load_ms(make_pg_matrix(d, "b.txt"))
  csv <- load_ms(make_pg_matrix(d, "c.csv", sep = ","))

  expect_equal(txt$wide, tsv$wide)
  expect_equal(csv$wide, tsv$wide)
})

test_that("n_annotation_cols overrides the annotation-column auto-detection", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  p <- make_pg_matrix(d)
  ms <- load_ms(p, n_annotation_cols = 2)    # only Protein.Group and Protein.Names

  # Genes and the description now count as sample columns -> 7 "samples", mostly NA
  expect_equal(ncol(ms$wide), 7L)
  expect_equal(colnames(ms$protein_meta), c("Protein.Group", "Protein.Names"))
})

test_that("duplicate sample columns and duplicate protein groups are made unique", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  p <- file.path(d, "dup.tsv")
  df <- data.frame(Protein.Group = c("P1", "P1", "P2"),
                   Protein.Names = "n", Genes = "g", First.Protein.Description = "d",
                   check.names = FALSE, stringsAsFactors = FALSE)
  df[["7_Plate 1"]] <- c(2, 4, 8)
  df[["7_Plate 2"]] <- c(2, 4, 8)    # same SampleID on a second plate
  utils::write.table(df, p, sep = "\t", quote = FALSE, row.names = FALSE)
  ms <- load_ms(p)

  expect_equal(colnames(ms$wide), c("7", "7_1"))
  expect_equal(rownames(ms$wide), c("P1", "P1_dup1", "P2"))
  expect_false(anyDuplicated(rownames(ms$wide)) > 0)
})

test_that("columns without the _Plate N suffix keep their name and get PlateID NA", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  p <- file.path(d, "noplate.tsv")
  df <- data.frame(Protein.Group = c("P1", "P2"), Protein.Names = "n", Genes = "g",
                   First.Protein.Description = "d", stringsAsFactors = FALSE)
  df$sampleA <- c(2, 4); df$sampleB <- c(8, 16)
  utils::write.table(df, p, sep = "\t", quote = FALSE, row.names = FALSE)
  sm <- load_ms(p)$sample_meta

  expect_equal(sm$SampleID, c("sampleA", "sampleB"))
  expect_true(all(is.na(sm$PlateID)))
})

test_that("load_massspec rejects missing files, bad formats and unknown sources", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  p <- make_pg_matrix(d)

  expect_error(load_ms(file.path(d, "nope.tsv")), "not found")
  bad <- file.path(d, "x.parquet"); writeLines("x", bad)
  expect_error(load_ms(bad), "Unsupported format")
  expect_error(load_ms(p, metadata_file = file.path(d, "nope.csv")), "metadata_file not found")
  expect_error(load_ms(p, metadata_file = make_pg_metadata(d), sample_col = "Nope"),
               "sample_col")
  # source is validated before anything is read
  expect_error(ELEUTHIA_load_massspec(p, source = "maxquant", verbose = FALSE), "should be")
})

test_that("the default source is dia-nn and dispatching it explicitly gives the same object", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  p <- make_pg_matrix(d)
  a <- load_ms(p)
  b <- load_ms(p, source = "dia-nn")

  expect_identical(a, b)
})

test_that("print.massspec_data summarises the object", {
  d <- tempfile("ms_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  ms <- load_ms(make_pg_matrix(d))
  out <- utils::capture.output(res <- print(ms))

  expect_identical(res, ms)       # returns its input invisibly
  expect_true(any(grepl("Proteins : 6", out)))
  expect_true(any(grepl("Samples  : 5", out)))
  expect_true(any(grepl("SAMPLE=4", out)))
})

# ---- export to the sample sheet -----------------------------------------------

# A processed object: group lives in sample_meta before export.
make_ms_for_export <- function() {
  d <- tempfile("ms_src_"); dir.create(d)
  ms <- load_ms(make_pg_matrix(d), metadata_file = make_pg_metadata(d,
          ids = c("1001", "1002", "SEP001-3", "SEP002-1")))
  unlink(d, recursive = TRUE)
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
