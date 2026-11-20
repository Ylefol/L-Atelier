# Shared fixture: simulated time-series counts with a known answer.
#
# 150 genes, 2 groups (ctrl / trt) x 3 timepoints x 4 replicates = 24 samples.
#   genes 1-10   rise at T3 in BOTH groups                  (a temporal effect)
#   genes 11-20  rise at T2 and T3 in 'trt' only            (a group x time effect)
#   genes 21-150 unchanged
make_timeseries <- function(seed = 11, n_genes = 150) {
  set.seed(seed)
  group <- rep(c("ctrl", "trt"), each = 12)
  tp    <- rep(rep(c("T1", "T2", "T3"), each = 4), 2)
  n     <- length(group)
  depth <- runif(n, 0.7, 1.4)
  base  <- rgamma(n_genes, 2, 0.02) + 20

  fc <- matrix(1, n_genes, n)
  fc[1:10, tp == "T3"] <- 4
  fc[11:20, group == "trt" & tp != "T1"] <- 4

  mu <- outer(base, depth) * fc
  counts <- matrix(rnbinom(length(mu), mu = mu, size = 20), n_genes)
  genes <- sprintf("g%03d", seq_len(n_genes))
  dimnames(counts) <- list(genes, sprintf("s%02d", seq_len(n)))

  list(counts = counts,
       targets = data.frame(group = group, timepoint = tp, row.names = colnames(counts),
                            stringsAsFactors = FALSE),
       temporal = genes[1:10], group_time = genes[11:20], null = genes[21:n_genes])
}

# DESeq2 reports dispersion-fit notes through message(); keep test logs readable.
quiet_msg <- function(expr) suppressMessages(quiet(expr))

sig_genes <- function(res, alpha = 0.05) {
  r <- res$results
  r$feature_id[!is.na(r$padj) & r$padj < alpha]
}
