# 04_runtime_scaling.R
# ------------------------------------------------------------
# Description: 
#   Benchmarks RSpectra (naive) vs rARPACK (function-based) 
#   eigendecomposition runtimes for varying sample sizes. 
# ============================================================

suppressPackageStartupMessages({
  library(RSpectra)
  library(rARPACK)
})

# Simulate genotype matrix (0,1,2)
simulate_genotypes <- function(n, m) {
  matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)
}

# Compute A_min for scaling
compute_Amin <- function(X) {
  n <- nrow(X)
  m <- ncol(X)
  X1 <- sweep(X, 2, colMeans(X))
  A <- (1 / m) * tcrossprod(X1) - matrix(1, n, n) / n
  min(A)
}

# RSpectra version (explicit matrix)
run_rspectra <- function(X, k) {
  n <- nrow(X)
  Phi <- popkin::popkin(X)
  t0 <- Sys.time()
  eigs <- RSpectra::eigs_sym(Phi, k = min(k, n - 1))
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  runtime
}

# rARPACK version (function-based)
run_rarpack <- function(X, k) {
  n <- nrow(X); m <- ncol(X)
  A_min <- compute_Amin(X)
  ones_n <- rep(1, n)
  X1 <- sweep(X, 2, colMeans(X))
  
  Afun <- function(v, args = NULL) {
    v <- as.numeric(v)
    term1 <- ones_n * sum(v) / n
    term2 <- (1 / m) * (X1 %*% (t(X1) %*% v)) - ones_n * sum(v) / n
    term1 - (1 / A_min) * term2
  }
  
  t0 <- Sys.time()
  eigs <- rARPACK::eigs_sym(Afun, k = min(k, n - 1), n = n)
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  runtime
}

# Benchmark loop
set.seed(123)
sizes <- c(10, 100, 1000)
k <- 10
results <- data.frame(n = sizes, RSpectra = NA, rARPACK = NA)

for (i in seq_along(sizes)) {
  n <- sizes[i]
  m <- 100 * n
  cat("\nRunning benchmarks for n =", n, ", m =", m, "...\n")
  X <- simulate_genotypes(n, m)
  results$RSpectra[i] <- run_rspectra(X, k)
  results$rARPACK[i]  <- run_rarpack(X, k)
}

# Save results 
dir.create("output", showWarnings = FALSE)
write.csv(results, "output/runtime_scaling.csv", row.names = FALSE)
cat("\nRuntime results saved to 'output/runtime_scaling.csv'\n")

# Plot for visualization 
png("output/runtime_scaling.png", width = 700, height = 500)
plot(results$n, results$RSpectra, type = "b", pch = 19, log = "xy",
     col = "red", xlab = "Sample size (n, log scale)",
     ylab = "Runtime (seconds, log scale)",
     main = "Runtime Scaling: RSpectra vs rARPACK")
lines(results$n, results$rARPACK, type = "b", pch = 19, col = "blue")
legend("topleft", legend = c("RSpectra", "rARPACK"),
       col = c("red", "blue"), pch = 19, bty = "n")
dev.off()

cat("Plot saved to 'output/runtime_scaling.png'\n")





