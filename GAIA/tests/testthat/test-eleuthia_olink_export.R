# ELEUTHIA_export_olink_qc: writes the filtered data, the QC log tables, plots and
# a summary built from olink_data$qc_log (nothing is assembled by hand). Uses the
# simulated Olink data from helper-olink.R, so the removed IDs have known answers.

export_fixture <- function() {
  ol <- run_olink_chain(load_sim_olink())
  assoc <- quiet(ARTEMIS_pc_metadata_association(
    matrix = ol$wide, sample_meta = ol$sample_meta, sample_col = "SampleID",
    n_pcs = 3, verbose = FALSE))
  qc <- quiet(AETHER_plot_olink_qc(ol, plate_col = "PlateID",
                                   plot_warn_proteins = TRUE, title = "t"))
  list(ol = ol, assoc = assoc, qc = qc)
}

do_export <- function(fx, out, ...) {
  quiet(ELEUTHIA_export_olink_qc(out, fx$ol, assoc = fx$assoc,
                                 qc_plots = list(pre = fx$qc, post = fx$qc),
                                 top_n_vars = 2, ...))
}

test_that("export writes the data, QC log tables, plots and summary, and returns the file list", {
  fx <- export_fixture()
  out <- tempfile("olink_out_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  files <- do_export(fx, out)

  expect_type(files, "character")
  expect_true(all(file.exists(files)))
  rel <- sub(paste0(out, "/"), "", files, fixed = TRUE)
  expect_true(all(c("data/olink_filtered.rds", "data/olink_sample_meta.csv",
                    "data/olink_assay_meta.csv", "qc/olink_qc_log.csv",
                    "qc/olink_removed_samples.csv", "qc/olink_removed_proteins.csv",
                    "qc/olink_flagged_samples.csv", "qc/olink_summary.txt",
                    "qc/olink_pre_npx_distributions.png", "qc/olink_post_sample_qc.csv",
                    "association/olink_pc_metadata_heatmap.png") %in% rel))
  expect_true(any(grepl("^pca/olink_pca_.*\\.png$", rel)))
})

test_that("the filtered RDS keeps the QC log, so it can be read back later", {
  fx <- export_fixture()
  out <- tempfile("olink_out_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  do_export(fx, out)

  back <- readRDS(file.path(out, "data", "olink_filtered.rds"))
  expect_s3_class(back, "olink_data")
  expect_equal(HADES_get_olink_qc_log(back)$step, HADES_get_olink_qc_log(fx$ol)$step)
})

test_that("the QC log CSVs match the log held in the object", {
  fx <- export_fixture()
  out <- tempfile("olink_out_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  do_export(fx, out)

  steps <- utils::read.csv(file.path(out, "qc", "olink_qc_log.csv"), stringsAsFactors = FALSE)
  expect_equal(steps$step, HADES_get_olink_qc_log(fx$ol)$step)
  expect_equal(nrow(steps), 6L)

  rs <- utils::read.csv(file.path(out, "qc", "olink_removed_samples.csv"), stringsAsFactors = FALSE)
  expect_setequal(rs$sample_id, c("X1", "N1", "N2", "S21", "S22", "S05"))
  rp <- utils::read.csv(file.path(out, "qc", "olink_removed_proteins.csv"), stringsAsFactors = FALSE)
  expect_setequal(rp$OlinkID, c("OID001", "OID002", "OID003"))
  fl <- utils::read.csv(file.path(out, "qc", "olink_flagged_samples.csv"), stringsAsFactors = FALSE)
  expect_equal(fl$sample_id, "S05")
})

test_that("the summary text is built from the log: every step, its arguments and its removed IDs", {
  fx <- export_fixture()
  out <- tempfile("olink_out_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  do_export(fx, out)
  txt <- readLines(file.path(out, "qc", "olink_summary.txt"))

  expect_true(any(grepl("FILTER CHAIN", txt)))
  for (s in c("HADES_filter_olink", "HADES_filter_olink_proteins",
              "HADES_filter_olink_assay_warn", "HADES_detect_outliers_olink"))
    expect_true(any(grepl(s, txt, fixed = TRUE)), info = s)
  expect_true(any(grepl("max_na_fraction=0.2", txt, fixed = TRUE)))
  expect_true(any(grepl("flagged", txt)))
  expect_true(any(grepl("S21", txt)))               # a removed sample ID
  expect_true(any(grepl("GENE001", txt)))           # a removed protein, by assay name
  expect_true(any(grepl("Proteins : 27", txt)))     # final dimensions
  expect_true(any(grepl("Samples  : 19", txt)))
})

test_that("prefix, plot_format, save_plots and save_rds control what is written", {
  fx <- export_fixture()

  out <- tempfile("olink_out_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  # The default pdf device cannot draw the em-dash / >= characters in the plot
  # titles (set by the AETHER plots) and warns; that is not what is tested here.
  files <- suppressWarnings(do_export(fx, out, prefix = "study", plot_format = "pdf"))
  rel <- sub(paste0(out, "/"), "", files, fixed = TRUE)
  expect_true("data/study_filtered.rds" %in% rel)
  expect_true("association/study_pc_metadata_heatmap.pdf" %in% rel)
  expect_false(any(grepl("olink_", basename(rel))))

  out2 <- tempfile("olink_out2_"); on.exit(unlink(out2, recursive = TRUE), add = TRUE)
  files2 <- do_export(fx, out2, save_plots = FALSE, save_rds = FALSE)
  rel2 <- sub(paste0(out2, "/"), "", files2, fixed = TRUE)
  expect_false(any(grepl("\\.(png|pdf)$", rel2)))
  expect_false(any(grepl("\\.rds$", rel2)))
  expect_true("qc/olink_qc_log.csv" %in% rel2)       # tables are still written
  expect_true("qc/olink_summary.txt" %in% rel2)
})

test_that("the association part is optional", {
  fx <- export_fixture()
  out <- tempfile("olink_out_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  files <- quiet(ELEUTHIA_export_olink_qc(out, fx$ol))

  expect_false(dir.exists(file.path(out, "association")))
  expect_false(dir.exists(file.path(out, "pca")))
  expect_true(file.exists(file.path(out, "qc", "olink_summary.txt")))
  txt <- readLines(file.path(out, "qc", "olink_summary.txt"))
  expect_false(any(grepl("PC-METADATA ASSOCIATION", txt)))
})

test_that("batch plots are saved under their names", {
  fx <- export_fixture()
  out <- tempfile("olink_out_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
  files <- quiet(ELEUTHIA_export_olink_qc(out, fx$ol, batch_plots = list(PlateID = p)))

  expect_true(file.exists(file.path(out, "qc", "olink_pca_batch_PlateID.png")))
})

test_that("an object without a QC log still exports, and the summary says so", {
  fx <- export_fixture()
  fx$ol$qc_log <- NULL
  out <- tempfile("olink_out_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  files <- quiet(ELEUTHIA_export_olink_qc(out, fx$ol))

  expect_false(file.exists(file.path(out, "qc", "olink_qc_log.csv")))
  txt <- readLines(file.path(out, "qc", "olink_summary.txt"))
  expect_true(any(grepl("No QC log is recorded", txt)))
  expect_true(any(grepl("Proteins : 27", txt)))
})

test_that("export validates its inputs", {
  fx <- export_fixture()
  out <- tempfile("olink_out_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)

  expect_error(ELEUTHIA_export_olink_qc(out, list()), "olink_data")
  expect_error(ELEUTHIA_export_olink_qc(out, fx$ol, assoc = list()), "artemis_pc_assoc")
  expect_error(ELEUTHIA_export_olink_qc(out, fx$ol, qc_plots = "x"), "qc_plots")
  expect_error(ELEUTHIA_export_olink_qc(out, fx$ol, batch_plots = list(1)), "named list")
})

test_that("the old hand-built qc_steps argument is gone", {
  fx <- export_fixture()
  out <- tempfile("olink_out_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)

  expect_false("qc_steps" %in% names(formals(ELEUTHIA_export_olink_qc)))
  expect_error(ELEUTHIA_export_olink_qc(out, fx$ol, qc_steps = list()), "unused argument")
})
