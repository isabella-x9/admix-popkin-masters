# 04_runtime_scaling.R
# ------------------------------------------------------------
# Description: 
#   Benchmarks RSpectra (naive) vs rARPACK (function-based) 
#   eigendecomposition runtimes for varying sample sizes. 
# ============================================================

suppressPackageStartupMessages({
  library(RSpectra)
  library(rARPACK) 
  library(popkin)
})

# Read n from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Please provide n")
n <- as.integer(args[1])
cat("Running for n =", n, "\n") 

# Fixed m 
m <- 10000
k <- 10

# Simulate genotype matrix (0,1,2)
simulate_genotypes <- function(n, m) {
  matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)
}

# RSpectra version (explicit matrix)
run_rspectra <- function(X, k) {
  t0 <- Sys.time()
  Phi <- popkin::popkin(X, loci_on_cols = TRUE)
  Theta <- inbr_diag(Phi)
  eigs <- RSpectra::eigs_sym(Theta, k = min(k, nrow(X) - 1))
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
}

# rARPACK version (function-based)
run_rarpack <- function(X, k) {
  n <- nrow(X)
  A_min <- min(popkin_A(X, loci_on_cols = TRUE)$A)
  X1 <- X - 1   # subtract 1 from every entry
  d <- rowMeans(X1^2 - 1) / A_min
  args <- list(X1 = X1, A_min = A_min, d = d) 
  
  source("Phi_prod.R")
  
  t0 <- Sys.time()
  eigs <- rARPACK::eigs_sym(Theta_prod, k = min(k, n - 1), n = n, args = args)
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  runtime
}


cat("\nRunning benchmarks for n =", n, ", m =", m, "...\n")
X <- simulate_genotypes(n, m)

time_RS <- run_rspectra(X, k)
time_RA <- run_rarpack(X, k)

# Save CSV for this n
dir.create("output", showWarnings = FALSE)
outfile <- sprintf("output/runtime_scaling_n_%d.csv", n)
write.csv(
  data.frame(n = n, m = m, RSpectra = time_RS, rARPACK = time_RA),
  outfile,
  row.names = FALSE
)
cat("Saved:", outfile, "\n")

# Plot for visualization 
pdf(sprintf("output/runtime_scaling_n_%d.pdf", n), width = 7, height = 5)
plot(c(n), c(time_RA),
     log = "xy", pch = 19, col = "red",
     xlab = "Sample size (n, log scale)",
     ylab = "Runtime (seconds, log scale)",
     main = sprintf("Runtime Scaling for n = %d", n),
     xlim = c(n/2, n*2), ylim = c(min(time_RA, time_RS)/2, max(time_RA, time_RS)*2)
)
points(c(n), c(time_RS), pch = 19, col = "blue")
legend("topleft", legend = c("rARPACK", "RSpectra"),
       col = c("red", "blue"), pch = 19)
dev.off()

cat("Saved:", sprintf("output/runtime_scaling_n_%d.pdf", n), "\n")
cat("\nDone.\n")





