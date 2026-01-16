suppressPackageStartupMessages({
  library(popkin)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/04b_compute_phi_amin.R <n>")

n <- as.integer(args[1])

geno_path <- sprintf("output/tmp/geno_n_%d.tsv", n)
phi_path  <- sprintf("output/tmp/Phi_n_%d.tsv", n)
amin_path <- sprintf("output/tmp/Phi_n_%d_Amin.rds", n)

stopifnot(file.exists(geno_path))
dir.create("output/tmp", showWarnings = FALSE, recursive = TRUE)

if (file.exists(phi_path) && file.exists(amin_path)) {
  cat("Exists, skipping Phi/Amin for n =", n, "\n")
  quit(save="no", status=0)
}

X <- as.matrix(read.table(geno_path, header=TRUE))
Phi <- popkin::popkin(X, loci_on_cols = TRUE)
A_min <- min(popkin_A(X, loci_on_cols = TRUE)$A)

write.table(Phi, phi_path, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
saveRDS(A_min, amin_path)

cat("Wrote:", phi_path, "\n")
cat("Wrote:", amin_path, "\n")