# 03_eigs_rarpack.R
# ------------------------------------------------------------
# Description: 
#   Scalable eigendecomposition of the Popkin kinship matrix using 
#   rARPACK, without forming Phi explicitly. 
#   It defines a matrix–vector product function based on the 
#   Popkin-EVD formulation and benchmarks runtime versus RSpectra. 
# ============================================================

suppressPackageStartupMessages({
  library(popkin)
  library(genio)
  library(rARPACK)
  library(argparser)
})

# Parse command-line arguments
p <- arg_parser("rARPACK eigendecomposition benchmark")
p <- add_argument(p, "--geno", help = "Path to genotype matrix", default = "data/test_geno.tsv")
p <- add_argument(p, "--amin", help = "Path to A_min RDS file",   default = "output/Phi_Amin.rds")
p <- add_argument(p, "--k",    help = "Number of top eigenvectors", default = 5)
argv <- parse_args(p)

# Load genotype matrix
cat("Loading genotype matrix from:", argv$geno, "\n")
X <- as.matrix(read.table(argv$geno, header = TRUE))
n <- nrow(X)
m <- ncol(X)
cat("Matrix dimensions:", n, "x", m, "\n") 

# Center columns (SNP means)
X1 <- X - 1   # subtract 1 from every entry

# Load A_min
cat("Loading A_min from:", argv$amin, "\n")
A_min <- readRDS(argv$amin)
cat("A_min =", round(A_min, 6), "\n")

# Define d 
d <- rowMeans( X1^2 - 1 ) / A_min 

# The Phi * v function already defined
source("scripts/Phi_prod.R")

# Compute top-k eigenvectors using rARPACK
k <- 5
start_time <- Sys.time()
end_time <- Sys.time()
args <- list(
  X1 = X1, A_min = A_min, d = d 
)
eigs <- eigs_sym(Theta_prod, k = k, n = n, args = args)
runtime <- round(difftime(Sys.time(), start_time, units = "secs"), 3)
cat("Eigendecomposition completed in", runtime, "seconds.\n")

# Save results
dir.create("output", showWarnings = FALSE)
eig_vals <- data.frame(eigenvalue = eigs$values)
eig_vecs <- as.data.frame(eigs$vectors)
colnames(eig_vecs) <- paste0("PC", seq_len(ncol(eig_vecs)))


write.table(eig_vals, "output/rarpack_vals.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(eig_vecs, "output/rarpack_vecs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(runtime, "output/rarpack_runtime.rds")
cat("Saved eigenvalues, eigenvectors, and runtime to 'output/'\n")





