# KHALKOS: shared helpers used by every other CAULDRON module.

test_that("validate_sce accepts a valid object and reports every missing piece", {
  sce <- make_test_sce()

  expect_true(KHALKOS_validate_sce(sce, require_assays = "counts",
                                   require_coldata = c("cell_type", "sample_id"),
                                   verbose = FALSE))

  err <- expect_error(
    KHALKOS_validate_sce(sce, require_assays = c("counts", "logcounts"),
                         require_reduceddims = "PCA",
                         require_coldata = "nope", require_metadata = "hvg",
                         verbose = FALSE),
    "SCE validation failed"
  )
  # all four problems are listed in one message, not just the first
  expect_match(conditionMessage(err), "logcounts")
  expect_match(conditionMessage(err), "PCA")
  expect_match(conditionMessage(err), "nope")
  expect_match(conditionMessage(err), "hvg")
})

test_that("validate_sce rejects objects that are not a SingleCellExperiment", {
  expect_error(KHALKOS_validate_sce(matrix(1:4, 2)), "not a SingleCellExperiment")
})

test_that("get_assay / set_assay round-trip and fail clearly", {
  sce <- make_test_sce()
  m <- KHALKOS_get_assay(sce, "counts")
  expect_equal(dim(m), dim(sce))

  sce2 <- KHALKOS_set_assay(sce, "doubled", m * 2)
  expect_true("doubled" %in% SummarizedExperiment::assayNames(sce2))
  expect_equal(KHALKOS_get_assay(sce2, "doubled"), m * 2)

  expect_error(KHALKOS_get_assay(sce, "missing"), "not found.*Available: counts")
  expect_null(KHALKOS_get_assay(sce, "missing", error_if_missing = FALSE))
})

test_that("default_palette returns the requested number of colours", {
  for (type in c("categorical", "sequential", "diverging")) {
    cols <- KHALKOS_default_palette(7, type)
    expect_length(cols, 7)
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", cols)))
  }
  # more categories than the curated 20: interpolated, not recycled or truncated
  expect_length(KHALKOS_default_palette(35), 35)
  expect_error(KHALKOS_default_palette(3, "nonsense"))
})

test_that("the categorical palette is stable: colour i never depends on n", {
  # Plots rely on this so a category keeps its colour whatever the subset size.
  expect_identical(KHALKOS_default_palette(5), KHALKOS_default_palette(12)[1:5])
})

test_that("assign_metadata_colours maps categorical columns only", {
  sce <- make_test_sce()
  sce$score <- seq_len(ncol(sce))   # numeric: must be skipped
  cols <- KHALKOS_assign_metadata_colours(sce, verbose = FALSE)

  expect_setequal(names(cols), c("cell_type", "sample_id"))
  expect_equal(names(cols$cell_type), c("A", "B", "C"))
  expect_equal(names(cols$sample_id), c("s1", "s2"))
  expect_false(anyDuplicated(unname(cols$cell_type)) > 0)
})

test_that("assign_metadata_colours is deterministic and sorts numeric-looking levels numerically", {
  sce <- make_test_sce()
  sce$cluster <- rep(c("10", "2", "1"), length.out = ncol(sce))
  a <- KHALKOS_assign_metadata_colours(sce, verbose = FALSE)
  b <- KHALKOS_assign_metadata_colours(sce, verbose = FALSE)

  expect_identical(a, b)
  expect_equal(names(a$cluster), c("1", "2", "10"))   # not "1", "10", "2"
})

test_that("assign_metadata_colours honours overrides, palettes and max_levels", {
  sce <- make_test_sce()

  ov <- KHALKOS_assign_metadata_colours(sce, overrides = list(cell_type = c(B = "#000000")),
                                        verbose = FALSE)
  expect_equal(unname(ov$cell_type["B"]), "#000000")

  pal <- KHALKOS_assign_metadata_colours(sce, palette = c("red", "blue", "green"),
                                         verbose = FALSE)
  expect_equal(unname(pal$cell_type), c("red", "blue", "green"))

  # sample_id has 2 levels, cell_type has 3: max_levels = 2 drops cell_type
  few <- KHALKOS_assign_metadata_colours(sce, max_levels = 2, verbose = FALSE)
  expect_equal(names(few), "sample_id")

  expect_error(KHALKOS_assign_metadata_colours(list()), "SingleCellExperiment")
})

test_that("check_packages errors, warns or passes depending on availability", {
  expect_true(KHALKOS_check_packages(c("stats", "utils")))

  expect_error(KHALKOS_check_packages("notARealPackage12345"),
               "not available: notARealPackage12345.*BiocManager::install")
  expect_error(KHALKOS_check_packages("notARealPackage12345",
                                      install_hint = "do this instead"),
               "do this instead")
  expect_warning(
    res <- KHALKOS_check_packages("notARealPackage12345", error = FALSE),
    "not available"
  )
  expect_false(res)
})

test_that("log_message routes by level", {
  expect_message(KHALKOS_log_message("hello"), "hello")
  expect_silent(KHALKOS_log_message("hello", verbose = FALSE))
  expect_warning(KHALKOS_log_message("careful", level = "warning"), "careful")
  expect_error(KHALKOS_log_message("broken", level = "error"), "broken")
})

test_that("theme_cauldron returns a ggplot2 theme usable on a plot", {
  th <- KHALKOS_theme_cauldron(base_size = 14)
  expect_s3_class(th, "theme")
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
    ggplot2::geom_point() + th
  expect_s3_class(p, "ggplot")
})
