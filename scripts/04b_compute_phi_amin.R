# scripts/04b_compute_phi_amin.R
suppressPackageStartupMessages({
  library(popkin)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/04b_compute_phi_amin.R <n>")

n <- as.integer(args[1])

geno_path <- sprintf("output/tmp/geno_n_%d.rds", n)
phi_path <- sprintf("output/tmp/Phi_n_%d.rds", n)
amin_path <- sprintf("output/tmp/Phi_n_%d_Amin.rds", n)
time_path <- sprintf("output/tmp/popkin_n_%d_timing.csv", n)

stopifnot(file.exists(geno_path))
dir.create("output/tmp", showWarnings = FALSE, recursive = TRUE)

# Skip only if ALL expected outputs exist (including timing)
if (file.exists(phi_path) && file.exists(amin_path) && file.exists(time_path)) {
  cat("Exists, skipping Phi/Amin/timing for n =", n, "\n")
  quit(save = "no", status = 0)
}

t_total0 <- Sys.time()

# Load X
t0 <- Sys.time()
X <- readRDS(geno_path) 
t_loadX <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# popkin()
t0 <- Sys.time()
Phi <- popkin::popkin(X, loci_on_cols = TRUE)
t_popkin <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# A_min via popkin_A()
t0 <- Sys.time()
A_min <- min(popkin::popkin_A(X, loci_on_cols = TRUE)$A) 
t_amin <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# Write outputs
t0 <- Sys.time()
saveRDS(Phi, phi_path) 
saveRDS(A_min, amin_path)
t_write <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

t_total <- as.numeric(difftime(Sys.time(), t_total0, units = "secs"))

# Save timing
write.csv(
  data.frame(
    n = n,
    load_X_sec = t_loadX,
    popkin_sec = t_popkin,
    popkin_A_sec = t_amin,
    write_outputs_sec = t_write,
    total_sec = t_total
  ),
  time_path,
  row.names = FALSE
)

cat("Wrote:", phi_path, "\n")
cat("Wrote:", amin_path, "\n")
cat("Wrote:", time_path, "\n")
