# 02_popkin_evd.R
# ------------------------------------------------------------
# Description: 
#   This script performs a scalable eigendecomposition of the Popkin kinship 
#   matrix (Phi) computed in Stage 1. 
#   It extracts the top-K eigenvalues and eigenvectors using an iterative 
#   Lanczos method in RSpectra. 
# ============================================================ 

# Load packages
suppressPackageStartupMessages({
  library(data.table)
  library(RSpectra)
})

# User inputs 
grm_path <- "output/Phi.tsv"     # path to Popkin GRM from Stage 1
k <- 5                          # number of top eigenvectors
out_prefix <- "output/eigen"

# Load Popkin kinship matrix
cat("Loading Popkin kinship matrix from:", grm_path, "\n")
Phi <- as.matrix(read.table(grm_path, header = TRUE, row.names = 1))
# Overwrite Phi by applying inbreeding diagonal (into coancestry matrix theta)
Phi <- inbr_diag(Phi) 
cat("Matrix dimensions:", paste(dim(Phi), collapse = " x "), "\n")

# Compute top-K eigenpairs
cat("Computing top", k, "eigenvectors using RSpectra::eigs_sym...\n")
start_time <- Sys.time()
eigs <- eigs_sym(Phi, k = k)
end_time <- Sys.time()

runtime <- round(difftime(end_time, start_time, units = "secs"), 3)
cat("Eigendecomposition completed in", runtime, "seconds.\n")

# Save results
dir.create("output", showWarnings = FALSE)
eig_vals <- data.frame(eigenvalue = eigs$values)
eig_vecs <- as.data.frame(eigs$vectors)
colnames(eig_vecs) <- paste0("PC", seq_len(k))

write.table(eig_vals, paste0(out_prefix, "_vals.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(eig_vecs, paste0(out_prefix, "_vecs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Eigenvalues and eigenvectors saved to 'output/'\n")



