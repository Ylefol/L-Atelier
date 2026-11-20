# Shared setup for the cross-package integration tests.
#
# These run AFTER install.R, against the installed packages (not the package
# sources), on every OS. They check that what HORIZON writes is what GAIA reads.
# Only simulated data is used; no external tools (samtools, STAR, ...) are needed.

# When run by hand from a development session, packages may already be attached
# (e.g. via pkgload::load_all); only attach them if they are not.
for (pkg in c("HORIZON", "GAIA")) {
  if (!paste0("package:", pkg) %in% search()) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

quiet <- function(expr) {
  # these functions report progress with cat(); keep the test log readable
  utils::capture.output(res <- force(expr))
  res
}

# Simulated bulk counts with a known answer: 10 genes up and 10 down in
# 'treated' (4x), 4 control + 4 treated samples with different sequencing
# depths, two batches.
simulate_counts <- function(n_genes = 200, n_per_group = 4, seed = 123) {
  set.seed(seed)
  n <- 2 * n_per_group
  condition <- rep(c("control", "treated"), each = n_per_group)
  depth <- rep(c(1, 1.5, 0.7, 2), length.out = n)

  base <- rgamma(n_genes, shape = 2, rate = 0.02) + 20
  fc <- rep(1, n_genes); fc[1:10] <- 4; fc[11:20] <- 1 / 4
  mu <- outer(base, depth)
  mu[, condition == "treated"] <- mu[, condition == "treated"] * fc

  counts <- matrix(rnbinom(length(mu), mu = mu, size = 20), nrow = n_genes)
  genes <- sprintf("gene%03d", seq_len(n_genes))
  samples <- sprintf("sample_%d", seq_len(n))
  dimnames(counts) <- list(genes, samples)

  list(counts = counts,
       condition = condition,
       batch = rep(c("b1", "b2"), length.out = n),
       up = genes[1:10], down = genes[11:20], null = genes[21:n_genes])
}

# Write per-sample count files exactly as HORIZON_run_count() does:
#   <output_dir>/<sample_id>/counts/<sample_id>_counts.txt
# tab-separated, header gene_id / count, integer counts. Returns a HORIZON
# sample sheet pointing at them.
write_horizon_run <- function(sim, output_dir, skip = character(0)) {
  for (sid in colnames(sim$counts)) {
    if (sid %in% skip) next
    d <- file.path(output_dir, sid, "counts")
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    df <- data.frame(gene_id = rownames(sim$counts),
                     count = as.integer(sim$counts[, sid]),
                     stringsAsFactors = FALSE)
    utils::write.table(df, file.path(d, paste0(sid, "_counts.txt")),
                       sep = "\t", quote = FALSE, row.names = FALSE)
  }
  data.frame(
    sample_id    = colnames(sim$counts),
    fastq_r1     = file.path(output_dir, paste0(colnames(sim$counts), "_R1.fastq.gz")),
    fastq_r2     = file.path(output_dir, paste0(colnames(sim$counts), "_R2.fastq.gz")),
    paired_end   = TRUE,
    output_dir   = output_dir,
    strandedness = "unstranded",
    condition    = sim$condition,
    batch        = sim$batch,
    stringsAsFactors = FALSE
  )
}

# The step a user performs by hand between the two packages: turn HORIZON's
# sample_metadata.csv into the 'targets' table ARTEMIS expects (sample IDs as
# row names).
targets_from_metadata <- function(metadata_csv) {
  meta <- utils::read.csv(metadata_csv, stringsAsFactors = FALSE)
  rownames(meta) <- meta$sample_id
  meta
}

run_dea <- function(counts, targets, group_col) {
  norm <- quiet(ARTEMIS_normalize_counts(counts, targets, group_col = group_col,
                                         verbose = FALSE))
  quiet(ARTEMIS_differential_counts(norm, reference = "control",
                                    experiment = "treated", verbose = FALSE))
}
