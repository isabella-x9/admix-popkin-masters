# scripts/04d_run_rarpack.R
# ------------------------------------------------------------
# rARPACK eigendecomposition for operator Theta
#
# Modes:
#   standard (default): complete simulated genotype matrix
#   miss: missingness-aware operator for real data
#
# Usage:
#   Rscript scripts/04d_run_rarpack.R <n>
#   Rscript scripts/04d_run_rarpack.R <n> standard
#   Rscript scripts/04d_run_rarpack.R <n> miss
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(rARPACK)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/04d_run_rarpack.R <n> [standard|miss]")
}

n <- as.integer(args[1])
mode <- if (length(args) >= 2) args[2] else "standard"

if (!mode %in% c("standard", "miss")) {
  stop("Mode must be 'standard' or 'miss'")
}

k <- min(10, n - 1)

geno_path <- sprintf("output/tmp/geno_n_%d.rds", n)
amin_path <- sprintf("output/tmp/Phi_n_%d_Amin.rds", n)

dir.create("output/tmp", showWarnings = FALSE, recursive = TRUE)

# Output prefixes depend on mode
out_prefix <- if (mode == "standard") {
  sprintf("output/tmp/rarpack_n_%d", n)
} else {
  sprintf("output/tmp/rarpack_miss_n_%d", n)
}

values_path  <- paste0(out_prefix, "_values.tsv")
vectors_path <- paste0(out_prefix, "_vectors.tsv")
rt_path      <- paste0(out_prefix, "_runtime.csv")
time_path    <- paste0(out_prefix, "_timing.csv")

# Skip if outputs already exist
if (file.exists(rt_path) && file.exists(time_path)) {
  cat("Outputs exist, skipping:", rt_path, "\n")
  quit(save = "no", status = 0)
}

stopifnot(file.exists(geno_path), file.exists(amin_path))

t_total0 <- Sys.time()

# ---- Load X ----
t0 <- Sys.time()
X <- readRDS(geno_path) 
t_loadX <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

if (nrow(X) != n) {
  stop("Loaded n = ", nrow(X), " but expected n = ", n)
}

# ---- Load A_min ----
t0 <- Sys.time()
A_min <- readRDS(amin_path)
t_loadAmin <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# ---- Source operator definitions ----
phi_prod_path <- file.path("scripts", "Phi_prod.R")
if (!file.exists(phi_prod_path)) stop("Missing: ", phi_prod_path)
source(phi_prod_path)

# ---- Build args depending on mode ----
if (mode == "standard") {
  
  # Complete-data centered matrix
  t0 <- Sys.time()
  X1 <- sweep(X, 2, colMeans(X))
  t_center <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  # Diagonal correction d
  t0 <- Sys.time()
  d <- rowMeans(X1^2 - 1) / A_min
  t_buildd <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  args_op <- list(X1 = X1, A_min = A_min, d = d)
  
  # ---- rARPACK eigs: standard operator ----
  t0 <- Sys.time()
  fit <- eigs_sym(Theta_prod, k = k, n = n, args = args_op)
  t_eigs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
} else {
  
  # Missingness-aware mode
  # Assumes Phi_prod.R defines Theta_prod_miss and that it knows
  # how to use the objects below
  
  # Simple centered matrix using column means with na.rm = TRUE
  t0 <- Sys.time()
  col_means <- colMeans(X, na.rm = TRUE)
  X1 <- sweep(X, 2, col_means)
  
  # Replace remaining NAs with 0 after centering so matrix products
  # can proceed
  X1[is.na(X1)] <- 0
  t_center <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  # Diagonal correction placeholder for missingness-aware mode
  t0 <- Sys.time()
  d <- rowMeans(X1^2 - 1) / A_min
  t_buildd <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  args_op <- list(X = X, X1 = X1, A_min = A_min, d = d)
  
  # ---- rARPACK eigs: missingness-aware operator ----
  t0 <- Sys.time()
  fit <- eigs_sym(Theta_prod_miss, k = k, n = n, args = args_op)
  t_eigs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
}

# ---- Write outputs ----
t0 <- Sys.time()
write.table(fit$values, values_path,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(fit$vectors, vectors_path,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
t_write <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

t_total <- as.numeric(difftime(Sys.time(), t_total0, units = "secs"))

# Runtime file: eigs only
write.csv(data.frame(runtime_sec = t_eigs), rt_path, row.names = FALSE)

# Detailed timing
write.csv(
  data.frame(
    n = n,
    k = k,
    mode = mode,
    load_X_sec = t_loadX,
    load_Amin_sec = t_loadAmin,
    center_X_sec = t_center,
    build_d_sec = t_buildd,
    eigs_sec = t_eigs,
    write_outputs_sec = t_write,
    total_sec = t_total
  ),
  time_path,
  row.names = FALSE
)

cat("Mode:", mode, "\n")
cat("rARPACK runtime (eigs-only, sec):", t_eigs, "\n")
cat("Wrote:", time_path, "\n")


