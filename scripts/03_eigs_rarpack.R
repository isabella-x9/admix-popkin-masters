# 03_runtime_benchmark.R
# ------------------------------------------------------------
# Description: 
#   Scalable eigendecomposition of the Pokin kinship matrix using 
#   rARPACK, without forming Phi explicitly 
#   Popkin kinship matrix function
# ============================================================

suppressPackageStartupMessages({
  library(popkin)
  library(genio)
  library(rARPACK)
})

# Load genotype matrix 
X <- as.matrix(read.table(argv$geno))
n <- nrow(X)
m <- ncol(X)

# Center columns (normalization part)
X1 <- X - rowMeans(X)

# Placeholder for 