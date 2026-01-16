# 04_runtime_scaling.R
# ------------------------------------------------------------
# Description:
#   Benchmark RSpectra vs rARPACK by calling the two separate scripts.
#   Uses runtime CSV files and produces one-row CSV per n.
# ============================================================

suppressPackageStartupMessages({
  library(popkin)
})

# Read n 
args <- commandArgs(trailingOnly = TRUE)

if (interactive()) {
  n <- 1000 
  cat("Interactive mode: using n =", n, "\n")
} else {
  if (length(args) == 0 || is.na(suppressWarnings(as.integer(args[1])))) {
    stop("Usage: Rscript scripts/04_runtime_scaling.R <n>")
  }
  n <- as.integer(args[1])
}

stopifnot(!is.na(n), n > 1)
cat("Running for n =", n, "\n")

# Fixed parameters
m <- 10000
k <- 10

# Output dirs
dir.create("output", showWarnings = FALSE)
dir.create("output/tmp", showWarnings = FALSE)

# Paths for this n
geno_path  <- sprintf("output/tmp/geno_n_%d.tsv", n)
phi_prefix <- sprintf("output/tmp/Phi_n_%d", n)
phi_path   <- paste0(phi_prefix, ".tsv")
amin_path  <- paste0(phi_prefix, "_Amin.rds")

rs_out <- sprintf("output/tmp/rspectra_n_%d", n)
ra_out <- sprintf("output/tmp/rarpack_n_%d", n)

# Simulate genotype matrix (0,1,2)
simulate_genotypes <- function(n, m) {
  matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)
}

cat("Simulating genotype matrix X (", n, "x", m, ")...\n", sep = "")
X <- simulate_genotypes(n, m)

# Save genotype to TSV for rARPACK script
write.table(X, geno_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cat("Saved genotype to:", geno_path, "\n")

# Compute Phi + A_min once (shared)
cat("Computing Phi + A_min once (shared inputs)...\n")
Phi <- popkin::popkin(X, loci_on_cols = TRUE)
A_min <- min(popkin_A(X, loci_on_cols = TRUE)$A)

write.table(Phi, phi_path, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
saveRDS(A_min, amin_path)

cat("Saved Phi to:", phi_path, "\n")
cat("Saved A_min to:", amin_path, "\n")

# Call RSpectra script 
cmd_rs <- sprintf(
  "Rscript scripts/03_eigs_rspectra.R --phi %s --k %d --out %s",
  shQuote(phi_path), k, shQuote(rs_out)
)
cat("\nRunning RSpectra:\n", cmd_rs, "\n", sep = "")
status_rs <- system(cmd_rs)
if (status_rs != 0) warning("RSpectra run failed for n = ", n)

# Call rARPACK script 
cmd_ra <- sprintf(
  "Rscript scripts/03_eigs_rarpack.R --geno %s --amin %s --k %d --out %s",
  shQuote(geno_path), shQuote(amin_path), k, shQuote(ra_out)
)
cat("\nRunning rARPACK:\n", cmd_ra, "\n", sep = "")
status_ra <- system(cmd_ra)
if (status_ra != 0) warning("rARPACK run failed for n = ", n)

# Read runtimes from CSV 
time_RS <- NA_real_
time_RA <- NA_real_

rs_rt_path <- paste0(rs_out, "_runtime.csv")
ra_rt_path <- paste0(ra_out, "_runtime.csv")

if (file.exists(rs_rt_path)) {
  time_RS <- read.csv(rs_rt_path)$runtime_sec[1]
}
if (file.exists(ra_rt_path)) {
  time_RA <- read.csv(ra_rt_path)$runtime_sec[1]
}

# Save per-n benchmark CSV
outfile <- sprintf("output/runtime_scaling_n_%d.csv", n)
write.csv(
  data.frame(n = n, m = m, k = k, RSpectra = time_RS, rARPACK = time_RA),
  outfile,
  row.names = FALSE
)
cat("\nSaved:", outfile, "\n")
cat("Done.\n")

