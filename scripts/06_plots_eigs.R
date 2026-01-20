# scripts/06_plot_eigs.R
# ------------------------------------------------------------
# Makes:
#   - line plot of eigenvalues (top k)
#   - line plot of eigenvectors (entries across individuals), colored by rank
#
# Uses outputs already written by:
#   output/tmp/rspectra_n_<n>_values.tsv, _vectors.tsv
#   output/tmp/rarpack_n_<n>_values.tsv, _vectors.tsv
#
# Notes:
#   - eigenvector signs are arbitrary; for display we "sign-align" each vector
#     by forcing its first entry to be positive.
# ------------------------------------------------------------

N_VALUES <- c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000)

pick_n <- function() {
  # prefer largest n where BOTH methods exist, else fall back to largest rARPACK
  both <- rev(N_VALUES)
  for (n in both) {
    rs_v <- sprintf("output/tmp/rspectra_n_%d_values.tsv", n)
    rs_u <- sprintf("output/tmp/rspectra_n_%d_vectors.tsv", n)
    ra_v <- sprintf("output/tmp/rarpack_n_%d_values.tsv", n)
    ra_u <- sprintf("output/tmp/rarpack_n_%d_vectors.tsv", n)
    if (file.exists(rs_v) && file.exists(rs_u) && file.exists(ra_v) && file.exists(ra_u)) return(n)
  }
  for (n in both) {
    ra_v <- sprintf("output/tmp/rarpack_n_%d_values.tsv", n)
    ra_u <- sprintf("output/tmp/rarpack_n_%d_vectors.tsv", n)
    if (file.exists(ra_v) && file.exists(ra_u)) return(n)
  }
  stop("No eig outputs found in output/tmp/. Run the pipeline first.")
}

read_vec <- function(path) as.numeric(scan(path, quiet = TRUE))
read_mat <- function(path) as.matrix(read.table(path, header = FALSE))

sign_align <- function(U) {
  # Force first entry positive for each column
  for (j in seq_len(ncol(U))) {
    if (U[1, j] < 0) U[, j] <- -U[, j]
  }
  U
}

n <- pick_n()
dir.create("output", showWarnings = FALSE, recursive = TRUE)

# rARPACK files (always expected)
ra_vals <- read_vec(sprintf("output/tmp/rarpack_n_%d_values.tsv", n))
ra_vecs <- sign_align(read_mat(sprintf("output/tmp/rarpack_n_%d_vectors.tsv", n)))

# RSpectra files (may not exist at n=100000)
rs_vals_path <- sprintf("output/tmp/rspectra_n_%d_values.tsv", n)
rs_vecs_path <- sprintf("output/tmp/rspectra_n_%d_vectors.tsv", n)
have_rs <- file.exists(rs_vals_path) && file.exists(rs_vecs_path)

if (have_rs) {
  rs_vals <- read_vec(rs_vals_path)
  rs_vecs <- sign_align(read_mat(rs_vecs_path))
}

k_ra <- length(ra_vals)
xk <- seq_len(k_ra)

# -------- eigenvalues line plot --------
pdf(sprintf("output/eigenvalues_n_%d.pdf", n), width = 7, height = 5)

plot(xk, ra_vals,
     type = "b", pch = 19,
     xlab = "Rank (1 = largest)",
     ylab = "Eigenvalue",
     main = sprintf("Top eigenvalues (n = %d)", n))

if (have_rs) {
  lines(seq_len(length(rs_vals)), rs_vals, type = "b", pch = 19)
  legend("topright", legend = c("rARPACK", "RSpectra"), pch = 19, lwd = 2)
} else {
  legend("topright", legend = c("rARPACK"), pch = 19, lwd = 2)
}

dev.off()

# -------- eigenvectors line plot (colored by rank) --------
# We plot entries across individuals (index 1...n)
# Negative values are allowed and meaningful
pdf(sprintf("output/eigenvectors_n_%d.pdf", n), width = 7.5, height = 5)

matplot(ra_vecs, type = "l", lty = 1,
        xlab = "Individual index",
        ylab = "Eigenvector entry",
        main = sprintf("Top eigenvectors from rARPACK (n = %d), color = rank", n))

legend("topright",
       legend = paste0("k=", seq_len(ncol(ra_vecs))),
       lty = 1, cex = 0.8)

dev.off()

# if RSpectra exists, also output its eigenvectors plot
if (have_rs) {
  pdf(sprintf("output/eigenvectors_rspectra_n_%d.pdf", n), width = 7.5, height = 5)
  
  matplot(rs_vecs, type = "l", lty = 1,
          xlab = "Individual index",
          ylab = "Eigenvector entry",
          main = sprintf("Top eigenvectors from RSpectra (n = %d), color = rank", n))
  
  legend("topright",
         legend = paste0("k=", seq_len(ncol(rs_vecs))),
         lty = 1, cex = 0.8)
  
  dev.off()
}

cat("Picked n =", n, "\n")
cat("Saved:\n",
    sprintf(" - output/eigenvalues_n_%d.pdf\n", n),
    sprintf(" - output/eigenvectors_n_%d.pdf\n", n), sep = "")
if (have_rs) cat(sprintf(" - output/eigenvectors_rspectra_n_%d.pdf\n", n))
