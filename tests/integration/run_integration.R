# Run the cross-package integration tests against the INSTALLED packages.
#   Rscript tests/integration/run_integration.R
# Fails (non-zero exit) if any test fails. Used by the install-check workflow.

if (!requireNamespace("testthat", quietly = TRUE)) install.packages("testthat")

cat("R:", R.version.string, "\nPlatform:", R.version$platform, "\n")
for (pkg in c("HORIZON", "GAIA")) {
  cat(pkg, "version:", as.character(utils::packageVersion(pkg)), "\n")
}

res <- testthat::test_dir("tests/integration", reporter = "summary",
                          stop_on_failure = TRUE, load_package = "none")
df <- as.data.frame(res)
cat(sprintf("\nIntegration tests: %d expectations, %d failed, %d errors, %d skipped, %d warnings\n",
            sum(df$nb), sum(df$failed), sum(df$error), sum(df$skipped), sum(df$warning)))
if (any(df$skipped)) stop("Some integration tests were skipped; none are expected to be.")
