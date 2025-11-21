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

# Simulate genotype matrix (0,1,2)
simulate_genotypes <- function(n, m) {
  matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)
}

# RSpectra version (explicit matrix)
run_rspectra <- function(X, k) {
  n <- nrow(X)
  t0 <- Sys.time()
  Phi <- popkin::popkin(X, loci_on_cols = TRUE)
  Theta <- inbr_diag(Phi) 
  eigs <- RSpectra::eigs_sym(Theta, k = min(k, n - 1))
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  runtime
} 

# rARPACK version (function-based)
run_rarpack <- function(X, k) {
  n <- nrow(X)
  A_min <- min(popkin_A(X, loci_on_cols = TRUE)$A)
  X1 <- X - 1   # subtract 1 from every entry
  d <- rowMeans( X1^2 - 1 ) / A_min 
  args <- list(
    X1 = X1, A_min = A_min, d = d 
  ) 
  
  source("scripts/Phi_prod.R")
  
  t0 <- Sys.time()
  eigs <- rARPACK::eigs_sym(Theta_prod, k = min(k, n - 1), n = n, args = args)
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  runtime
}

# Benchmark loop
set.seed(123)
sizes <- c(10, 100, 1000, 10000, 100000, 1000000)
k <- 10
results <- data.frame(n = sizes, RSpectra = NA, rARPACK = NA)

for (i in seq_along(sizes)) {
  n <- sizes[i]
  m <- 10000
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
plot(results$n, results$rARPACK, type = "b", pch = 19, log = "xy",
     col = "red", xlab = "Sample size (n, log scale)",
     ylab = "Runtime (seconds, log scale)",
     main = "Runtime Scaling: RSpectra vs rARPACK")
lines(results$n, results$RSpectra, type = "b", pch = 19, col = "blue")
legend("topleft", legend = c("rARPACK", "RSpectra"),
       col = c("red", "blue"), pch = 19, bty = "n")
dev.off()

cat("Plot saved to 'output/runtime_scaling.png'\n")





