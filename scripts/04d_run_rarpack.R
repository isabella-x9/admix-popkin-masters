# scripts/04d_run_rarpack.R
suppressPackageStartupMessages({
  library(rARPACK)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/04d_run_rarpack.R <n>")

n <- as.integer(args[1])
k <- min(10, n - 1)

geno_path <- sprintf("output/tmp/geno_n_%d.tsv", n)
amin_path <- sprintf("output/tmp/Phi_n_%d_Amin.rds", n)

out_prefix <- sprintf("output/tmp/rarpack_n_%d", n)
rt_path    <- paste0(out_prefix, "_runtime.csv")
time_path  <- paste0(out_prefix, "_timing.csv")

stopifnot(file.exists(geno_path), file.exists(amin_path))
dir.create("output/tmp", showWarnings = FALSE, recursive = TRUE)

# Skip if outputs exist
if (file.exists(rt_path) && file.exists(time_path)) {
  cat("rARPACK outputs exist, skipping:", rt_path, "\n")
  quit(save = "no", status = 0)
}

t_total0 <- Sys.time()

# Load X
t0 <- Sys.time()
X <- as.matrix(read.table(geno_path, header = TRUE))
t_loadX <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

n_check <- nrow(X)
if (n_check != n) stop("Loaded n = ", n_check, " but expected n = ", n)

p_obs <- rowMeans( !is.na( X ) )

# Center X -> X1
t0 <- Sys.time()
X1 <- X - 1
X1[ is.na( X1 ) ] <- 0
t_center <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# Load A_min
t0 <- Sys.time()
A_min <- readRDS(amin_path)
t_loadAmin <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# Build d
t0 <- Sys.time()
d <- rowMeans(X1^2 - 1) / A_min / p_obs
t_buildd <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# Phi_prod.R 
phi_prod_path <- file.path("scripts", "Phi_prod.R")
if (!file.exists(phi_prod_path)) stop("Missing: ", phi_prod_path)

source(phi_prod_path)

args_op <- list(X1 = X1, A_min = A_min, d = d)

# Eigs
t0 <- Sys.time()
fit <- eigs_sym(Theta_prod, k = k, n = n, args = args_op)
fit_miss <- eigs_sym(Theta_prod_miss, k = k, n = n, args = args_op)
t_eigs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# Write outputs
t0 <- Sys.time()
write.table(fit$values, paste0(out_prefix, "_values.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(fit$vectors, paste0(out_prefix, "_vectors.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
t_write <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

t_total <- as.numeric(difftime(Sys.time(), t_total0, units = "secs"))

# Save runtime 
write.csv(data.frame(runtime_sec = t_eigs), rt_path, row.names = FALSE)

# Save detailed timing
write.csv(
  data.frame(
    n = n,
    k = k,
    load_X_sec = t_loadX,
    center_X_sec = t_center,
    load_Amin_sec = t_loadAmin,
    build_d_sec = t_buildd,
    eigs_sec = t_eigs,
    write_outputs_sec = t_write,
    total_sec = t_total
  ),
  time_path,
  row.names = FALSE
)

cat("rARPACK runtime (sec):", t_eigs, "\n")
cat("rARPACK done.\n")
cat("Wrote:", time_path, "\n")


