suppressPackageStartupMessages({})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/04c_run_rspectra.R <n>")

n <- as.integer(args[1])
k <- 10

phi_path <- sprintf("output/tmp/Phi_n_%d.tsv", n)
out_prefix <- sprintf("output/tmp/rspectra_n_%d", n)
rt_path <- paste0(out_prefix, "_runtime.csv")

stopifnot(file.exists(phi_path))

if (file.exists(rt_path)) {
  cat("RSpectra runtime exists, skipping:", rt_path, "\n")
  quit(save="no", status=0)
}

cmd <- sprintf(
  "Rscript scripts/03_eigs_rspectra.R --phi %s --k %d --out %s",
  shQuote(phi_path), k, shQuote(out_prefix)
)

cat("Running:", cmd, "\n")
status <- system(cmd)
if (status != 0) stop("RSpectra failed for n = ", n) 

cat("RSpectra done.\n")
