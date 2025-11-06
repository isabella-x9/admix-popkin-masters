# 01_build_popkin_grm.R
# ------------------------------------------------------------
# Description:
#   This script calculates the Popkin kinship matrix (Phi) and marker weights (M)
#   from a general genotype matrix file, and stores A_min needed for later 
#   eigendecomposition stages 02 and 03. 
#   This is the first stage of the project. 
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
n <- nrow(X)
m <- ncol(X)
cat("Loaded genotype matrix of dimension", n, "x", m, "\n") 

# Compute Popkin kinship matrix
cat("Computing Popkin kinship matrix...\n")
Phi <- popkin(X) 
##### don't transpose it if it is a BEDMatrix object, because it is not a real matrix #########
cat("Kinship matrix dimensions:", paste(dim(Phi), collapse = " x "), "\n")

# Compute marker weights (mean allele frequency variance)
cat("Computing marker weights...\n")
p <- colMeans(X) / 2      # estimated allele frequencies
M <- 2 * p * (1 - p)      # variance per SNP 
cat("Mean marker variance:", mean(M), "\n")

# Compute A and A_min (for Popkin-EVD)
cat("Computing A and A_min for Popkin-EVD...\n")
X1 <- X - rowMeans(X)

# A = (1/m) * X1'X1 − 1_n 1_n'/n
ones_n <- matrix(1, n, n)
A <- (1 / m) * tcrossprod(X1) - ones_n / n
A_min <- min(A)
cat("Computed A_min:", A_min, "\n")

# Save outputs
dir.create("output", showWarnings = FALSE)
write.table(Phi, paste0(out_prefix, ".tsv"),
            sep = "\t", quote = FALSE)
write.table(M, paste0(out_prefix, "_marker_weights.tsv"),
            sep = "\t", quote = FALSE, col.names = FALSE)
saveRDS(A_min, paste0(out_prefix, "_Amin.rds"))

cat("Saved Phi, marker weights, and A_min to 'output/'\n")

cat("\nStage 01 complete.\n")






