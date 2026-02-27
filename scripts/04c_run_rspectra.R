# scripts/04c_run_rspectra.R
# ------------------------------------------------------------
# RSpectra eigendecomposition for explicit THETA (Theta)
#   Theta := Phi - diag(d)
# where d is computed the same way as in 04d_run_rarpack.R:
#   X1 <- sweep(X, 2, colMeans(X))
#   d  <- rowMeans(X1^2 - 1) / A_min
#
# Inputs:
#   output/tmp/Phi_n_<n>.tsv
#   output/tmp/geno_n_<n>.tsv
#   output/tmp/Phi_n_<n>_Amin.rds
#
# Outputs:
#   output/tmp/rspectra_n_<n>_values.tsv
#   output/tmp/rspectra_n_<n>_vectors.tsv
#   output/tmp/rspectra_n_<n>_runtime.csv
#   output/tmp/rspectra_n_<n>_timing.csv
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(RSpectra)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/04c_run_rspectra.R <n>")

n <- as.integer(args[1])
k <- min(10, n - 1)

phi_path   <- sprintf("output/tmp/Phi_n_%d.tsv", n)
geno_path  <- sprintf("output/tmp/geno_n_%d.tsv", n)
amin_path  <- sprintf("output/tmp/Phi_n_%d_Amin.rds", n)

out_prefix <- sprintf("output/tmp/rspectra_theta_n_%d", n)
rt_path    <- paste0(out_prefix, "_runtime.csv")
time_path  <- paste0(out_prefix, "_timing.csv")

dir.create("output/tmp", showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(phi_path), file.exists(geno_path), file.exists(amin_path))

# Skip if runtime exists
if (file.exists(rt_path) && file.exists(time_path)) {
  cat("RSpectra outputs exist, skipping:", rt_path, "\n")
  quit(save = "no", status = 0)
}

t_total0 <- Sys.time()

# ---- Load Phi ----
t0 <- Sys.time()
Phi <- as.matrix(read.table(phi_path, header = TRUE, row.names = 1))
t_loadPhi <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# ---- Load X ----
t0 <- Sys.time()
X <- as.matrix(read.table(geno_path, header = TRUE))
t_loadX <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

n_check <- nrow(X)
if (n_check != n) stop("Loaded n = ", n_check, " but expected n = ", n)

# ---- Load A_min ----
t0 <- Sys.time()
A_min <- readRDS(amin_path)
t_loadAmin <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# ---- Center X -> X1 ----
t0 <- Sys.time()
X1 <- sweep(X, 2, colMeans(X))
t_center <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# ---- Build d ----
t0 <- Sys.time()
d <- rowMeans(X1^2 - 1) / A_min
t_buildd <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# ---- Build Theta explicitly: Theta = Phi - diag(d) ----
t0 <- Sys.time()
Theta <- Phi
diag(Theta) <- diag(Theta) - d
t_buildTheta <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# ---- Eigs on Theta ----
t0 <- Sys.time()
fit <- eigs_sym(Theta, k = k, which = "LA")
t_eigs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# ---- Write outputs ----
t0 <- Sys.time()
write.table(fit$values, paste0(out_prefix, "_values.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(fit$vectors, paste0(out_prefix, "_vectors.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
t_write <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

t_total <- as.numeric(difftime(Sys.time(), t_total0, units = "secs"))

# Keep runtime as eigs only
write.csv(data.frame(runtime_sec = t_eigs), rt_path, row.names = FALSE)

# Save detailed timing
write.csv(
  data.frame(
    n = n,
    k = k,
    load_Phi_sec = t_loadPhi,
    load_X_sec = t_loadX,
    load_Amin_sec = t_loadAmin,
    center_X_sec = t_center,
    build_d_sec = t_buildd,
    build_Theta_sec = t_buildTheta,
    eigs_sec = t_eigs,
    write_outputs_sec = t_write,
    total_sec = t_total
  ),
  time_path,
  row.names = FALSE
)

cat("RSpectra runtime (eigs-only, sec):", t_eigs, "\n")
cat("RSpectra done.\n")
cat("Wrote:", time_path, "\n")
