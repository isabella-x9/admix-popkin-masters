# scripts/04e_combine_csv.R
# ------------------------------------------------------------
# Combine per-n totals into output/runtime_scaling_n_<n>.csv
# Uses timing CSVs so RSpectra and rARPACK are comparable (end-to-end).
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/04e_combine_csv.R <n>")

n <- as.integer(args[1])
m <- 10000
k <- 10

# Timing files 
popkin_time <- sprintf("output/tmp/popkin_n_%d_timing.csv", n)
rs_time     <- sprintf("output/tmp/rspectra_n_%d_timing.csv", n)
ra_time     <- sprintf("output/tmp/rarpack_n_%d_timing.csv", n)

# Fallback runtime files 
rs_rt <- sprintf("output/tmp/rspectra_n_%d_runtime.csv", n)
ra_rt <- sprintf("output/tmp/rarpack_n_%d_runtime.csv", n)

# Read totals
PhiTime <- if (file.exists(popkin_time)) read.csv(popkin_time)$total_sec[1] else NA_real_

RSpectraTotal <- if (file.exists(rs_time)) read.csv(rs_time)$total_sec[1] else {
  if (file.exists(rs_rt)) read.csv(rs_rt)$runtime_sec[1] else NA_real_
}

# If you want eigs-only too 
RSpectraEigs <- if (file.exists(rs_time)) read.csv(rs_time)$eigs_sec[1] else NA_real_

rARPACKTotal <- if (file.exists(ra_time)) read.csv(ra_time)$total_sec[1] else {
  if (file.exists(ra_rt)) read.csv(ra_rt)$runtime_sec[1] else NA_real_
}

dir.create("output", showWarnings = FALSE, recursive = TRUE)
outfile <- sprintf("output/runtime_scaling_n_%d.csv", n)

write.csv(
  data.frame(
    n = n, m = m, k = k,
    PhiTime = PhiTime,                 # bookkeeping
    RSpectraEigs = RSpectraEigs,       # bookkeeping
    RSpectraTotal = RSpectraTotal,  
    rARPACKTotal = rARPACKTotal   
  ),
  outfile,
  row.names = FALSE
)

cat("Wrote:", outfile, "\n")
