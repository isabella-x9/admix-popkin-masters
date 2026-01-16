suppressPackageStartupMessages({})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/04d_run_rarpack.R <n>")

n <- as.integer(args[1])
k <- 10

geno_path <- sprintf("output/tmp/geno_n_%d.tsv", n)
amin_path <- sprintf("output/tmp/Phi_n_%d_Amin.rds", n)
out_prefix <- sprintf("output/tmp/rarpack_n_%d", n)
rt_path <- paste0(out_prefix, "_runtime.csv")

stopifnot(file.exists(geno_path), file.exists(amin_path))

if (file.exists(rt_path)) {
  cat("Exists, skipping rARPACK for n =", n, "\n")
  quit(save="no", status=0)
}

cmd <- sprintf(
  "Rscript scripts/03_eigs_rarpack.R --geno %s --amin %s --k %d --out %s",
  shQuote(geno_path), shQuote(amin_path), k, shQuote(out_prefix)
)

cat("Running:", cmd, "\n")
status <- system(cmd)
if (status != 0) stop("rARPACK failed for n = ", n) 
