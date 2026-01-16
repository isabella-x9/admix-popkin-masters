suppressPackageStartupMessages({})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/04e_combine_csv.R <n>")

n <- as.integer(args[1])
m <- 10000
k <- 10

rs_rt <- sprintf("output/tmp/rspectra_n_%d_runtime.csv", n)
ra_rt <- sprintf("output/tmp/rarpack_n_%d_runtime.csv", n)

time_RS <- if (file.exists(rs_rt)) read.csv(rs_rt)$runtime_sec[1] else NA_real_
time_RA <- if (file.exists(ra_rt)) read.csv(ra_rt)$runtime_sec[1] else NA_real_

dir.create("output", showWarnings = FALSE, recursive = TRUE)
outfile <- sprintf("output/runtime_scaling_n_%d.csv", n)

write.csv(
  data.frame(n=n, m=m, k=k, RSpectra=time_RS, rARPACK=time_RA),
  outfile,
  row.names = FALSE
)

cat("Wrote:", outfile, "\n")