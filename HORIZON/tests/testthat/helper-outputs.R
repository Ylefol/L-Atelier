# Simulated per-sample HORIZON outputs, laid out the way HORIZON writes them:
#   <output_dir>/<sample_id>/counts/<sample_id>_counts.txt     (gene_id, count)
#   <output_dir>/<sample_id>/salmon/quant.sf                   (Name, ..., TPM, NumReads)
#   <output_dir>/<sample_id>/tecount/<sample_id>.cntTable      (id, count)
# No sequencing data or external tools are involved: only the files the
# aggregation functions read are written.

sim_sheet <- function(out, ids = c("s1", "s2", "s3")) {
  data.frame(sample_id = ids, output_dir = out,
             condition = rep(c("ctrl", "treat"), length.out = length(ids)),
             stringsAsFactors = FALSE)
}

# genes g1..g5 with a count per sample that encodes the sample index (s1 -> x1 ...)
write_sim_counts <- function(out, ids = c("s1", "s2", "s3"), genes = paste0("g", 1:5)) {
  for (i in seq_along(ids)) {
    d <- file.path(out, ids[i], "counts"); dir.create(d, recursive = TRUE, showWarnings = FALSE)
    utils::write.table(data.frame(gene_id = genes, count = i * 10L + seq_along(genes)),
                       file.path(d, paste0(ids[i], "_counts.txt")),
                       sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

write_sim_salmon <- function(out, ids = c("s1", "s2", "s3"), tx = paste0("tx", 1:4)) {
  for (i in seq_along(ids)) {
    d <- file.path(out, ids[i], "salmon"); dir.create(d, recursive = TRUE, showWarnings = FALSE)
    utils::write.table(data.frame(Name = tx, Length = 1000L, EffectiveLength = 900,
                                  TPM = i * 100 + seq_along(tx),
                                  NumReads = i * 1000 + seq_along(tx)),
                       file.path(d, "quant.sf"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

# A RepeatMasker-style TE GTF (gene_id / family_id / class_id attributes)
write_te_gtf <- function(path) {
  attrs <- c('gene_id "AluY"; family_id "Alu"; class_id "SINE";',
             'gene_id "L1PA2"; family_id "L1"; class_id "LINE";',
             'gene_id "AluY"; family_id "Alu"; class_id "SINE";')   # a repeated subfamily
  writeLines(paste("chr1", "rmsk", "exon", c(100, 500, 900), c(300, 800, 1200), ".", "+", ".", attrs,
                   sep = "\t"), path)
  path
}

# TEcount cntTable: gene rows and TE rows (composite id gene:family:class) mixed
write_sim_tecount <- function(out, ids = c("s1", "s2", "s3")) {
  for (i in seq_along(ids)) {
    d <- file.path(out, ids[i], "tecount"); dir.create(d, recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
      data.frame(`gene/TE` = c("ENSG1", "ENSG2", "AluY:Alu:SINE", "L1PA2:L1:LINE"),
                 count = c(5L, 6L, i * 10L, i * 20L), check.names = FALSE),
      file.path(d, paste0(ids[i], ".cntTable")), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

quiet <- function(expr) {
  utils::capture.output(res <- force(expr))
  res
}
