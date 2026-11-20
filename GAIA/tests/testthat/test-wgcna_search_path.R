# ARTEMIS_wgcna_detect_modules() needs WGCNA's own cor() on the search path
# (blockwiseModules otherwise fails with "unused arguments (weights.x = ...)").
# GAIA attaches WGCNA for the duration of the call only. These tests check that
# it works AND that the user's session is left exactly as it was found.

# Three planted co-expression modules (60 genes each) plus 60 noise genes.
simulate_modules <- function(n_samples = 60, genes_per_module = 60,
                             n_modules = 3, n_noise = 60, seed = 42) {
  set.seed(seed)
  block <- function() {
    f <- rnorm(n_samples)
    vapply(seq_len(genes_per_module),
           function(i) f * runif(1, 0.7, 1) + rnorm(n_samples, sd = 0.5),
           numeric(n_samples))
  }
  m <- do.call(cbind, c(replicate(n_modules, block(), simplify = FALSE),
                        list(matrix(rnorm(n_samples * n_noise), n_samples))))
  colnames(m) <- paste0("g", seq_len(ncol(m)))
  rownames(m) <- paste0("s", seq_len(n_samples))
  m
}

run_detect <- function(m) {
  ARTEMIS_wgcna_detect_modules(m, power = 6, min_module_size = 20,
                               verbose = FALSE)
}

test_that("module detection works when WGCNA is not attached, and cleans up", {
  skip_if(
    "package:WGCNA" %in% search(),
    "WGCNA is already attached in this session"
  )
  m <- simulate_modules()

  res <- suppressMessages(run_detect(m))

  # The planted structure is recovered: 3 modules, each made of one planted block
  expect_s3_class(res, "wgcna_modules")
  expect_equal(res$n_modules, 3)
  block_of <- rep(1:4, each = 60)   # 4 = noise genes
  for (b in 1:3) {
    mods <- unique(res$module_names[block_of == b])
    mods <- mods[mods != "module_0"]
    expect_length(mods, 1)
  }

  # The session is left as it was found
  expect_false("package:WGCNA" %in% search())
  expect_false(exists("cor", envir = globalenv(), inherits = FALSE))
  expect_identical(environmentName(environment(get("cor", globalenv()))), "stats")
})

test_that("an already attached WGCNA stays attached", {
  skip_if(
    "package:WGCNA" %in% search(),
    "WGCNA is already attached in this session"
  )
  suppressMessages(attachNamespace("WGCNA"))
  on.exit(detach("package:WGCNA", unload = FALSE), add = TRUE)

  res <- suppressMessages(run_detect(simulate_modules()))

  expect_equal(res$n_modules, 3)
  expect_true("package:WGCNA" %in% search())
})

test_that("search path is restored even when the analysis fails", {
  skip_if(
    "package:WGCNA" %in% search(),
    "WGCNA is already attached in this session"
  )
  m <- simulate_modules()
  # a non-numeric power makes the call fail after WGCNA has been attached
  expect_error(suppressWarnings(suppressMessages(
    ARTEMIS_wgcna_detect_modules(m, power = "not a number", verbose = FALSE)
  )), "Could not determine power")
  expect_false("package:WGCNA" %in% search())
})
