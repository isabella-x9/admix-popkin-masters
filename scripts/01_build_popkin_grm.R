# 01_build_popkin_grm.R
# ------------------------------------------------------------
# Description:
#   This script calculates the Popkin kinship matrix (Phi) and marker weights (M)
#   from a general genotype matrix file. This is the first stage of the project. 
# ============================================================

# Load packages
suppressPackageStartupMessages({
  library(popkin)
  library(genio)
})

# User inputs
geno_path <- "data/test_geno.tsv"   # path to genotype matrix
out_prefix <- "output/Phi"

# Load genotype matrix
cat("Loading genotype matrix from:", geno_path, "\n")
X <- as.matrix(read.table(geno_path, header = TRUE))
cat("Matrix dimensions:", paste(dim(X), collapse = " x "), "\n")

# Compute Popkin kinship matrix
cat("Computing Popkin kinship matrix...\n")
Phi <- popkin(t(X))
cat("Kinship matrix dimensions:", paste(dim(Phi), collapse = " x "), "\n")

# Compute marker weights (mean allele frequency variance)
cat("Computing marker weights...\n")
p <- colMeans(X) / 2                  # estimated allele frequencies
M <- 2 * p * (1 - p)                  # variance per SNP
cat("Mean marker variance:", mean(M), "\n")

# Save outputs
dir.create("output", showWarnings = FALSE)
write.table(Phi, paste0(out_prefix, ".tsv"), sep = "\t", quote = FALSE)
write.table(M, paste0(out_prefix, "_marker_weights.tsv"), sep = "\t", quote = FALSE, col.names = FALSE)
cat("Saved Popkin kinship matrix and marker weights to 'output/'\n")







