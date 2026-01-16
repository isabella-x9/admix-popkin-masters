suppressPackageStartupMessages({})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/04a_simulate_geno.R <n>")

n <- as.integer(args[1])
m <- 10000

dir.create("output/tmp", showWarnings = FALSE, recursive = TRUE)
geno_path <- sprintf("output/tmp/geno_n_%d.tsv", n)

if (file.exists(geno_path)) {
  cat("Exists, skipping:", geno_path, "\n")
  quit(save="no", status=0)
}

X <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)
colnames(X) <- paste0("SNP", seq_len(m))
write.table(X, geno_path, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

cat("Wrote:", geno_path, "\n")