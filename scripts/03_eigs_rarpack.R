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

# Load genotype matrix 
p <- arg_parser("rARPACK eigendecomposition benchmark")
p <- add_argument(p, "--geno", help = "Path to genotype matrix", default = "data/test_geno.tsv")
argv <- parse_args(p) 

X <- as.matrix(read.table(argv$geno, header = TRUE)) 
n <- nrow(X)
m <- ncol(X)

# Center columns (normalization part)
X1 <- sweep(X, 2, colMeans(X))  # center by SNP means 

# Load A_min from Stage 01
A_min <- readRDS("output/Phi_Amin.rds") 

# Define Phi * v function
ones_n <- rep(1, n)
Phi_prod <- function(v) {
  v <- as.numeric(v)
  term1 <- ones_n * sum(v) / n
  term2 <- (1 / m) * (X1 %*% (t(X1) %*% v)) - ones_n * sum(v) / n
  term1 - (1 / A_min) * term2
}

# Compute top-k eigenvectors using rARPACK
k <- 5
start_time <- Sys.time()
##########################################################################
eigs <- eigs_sym(Phi_prod, k = k, n = n)
########## Error in (function (v)  : unused argument (NULL) ###############

end_time <- Sys.time()
runtime <- round(difftime(end_time, start_time, units = "secs"), 3)
cat("Eigendecomposition completed in", runtime, "seconds.\n")

# Save results
eig_vals <- data.frame(eigenvalue = eigs$values)
eig_vecs <- as.data.frame(eigs$vectors)
colnames(eig_vecs) <- paste0("PC", seq_len(ncol(eig_vecs)))


write.table(eig_vals, "output/rarpack_vals.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(eig_vecs, "output/rarpack_vecs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(runtime, "output/rarpack_runtime.rds")
cat("Saved eigenvalues, eigenvectors, and runtime to 'output/'\n")





