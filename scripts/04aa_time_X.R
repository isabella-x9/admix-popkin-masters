# scripts/04aa_time_X.R
# ------------------------------------------------------------
# Time loading the genotype matrix X from RDS
#
# Input:
#   output/tmp/geno_n_<n>.rds
#
# Output:
#   output/tmp/X_n_<n>_timing.csv
# ------------------------------------------------------------

suppressPackageStartupMessages({})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/04aa_time_X.R <n>")

n <- as.integer(args[1])

geno_path <- sprintf("output/tmp/geno_n_%d.rds", n)
time_path <- sprintf("output/tmp/X_n_%d_timing.csv", n)

dir.create("output/tmp", showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(geno_path))

# Skip if timing file already exists
if (file.exists(time_path)) {
  cat("Exists, skipping:", time_path, "\n")
  quit(save = "no", status = 0)
}

t0 <- Sys.time()
X <- readRDS(geno_path)
t_load <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

write.csv(
  data.frame(
    n = n,
    m = ncol(X),
    load_X_sec = t_load
  ),
  time_path,
  row.names = FALSE
)

cat("Wrote:", time_path, "\n")