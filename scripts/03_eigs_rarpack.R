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

p <- arg_parser("rARPACK eigendecomposition (operator-based)")
p <- add_argument(p, "--geno", help="Path to genotype matrix (tsv)", default="data/test_geno.tsv")
p <- add_argument(p, "--amin", help="Path to A_min RDS file", default="output/Phi_Amin.rds")
p <- add_argument(p, "--k",    help="Number of top eigenvectors", type="integer", default=10)
p <- add_argument(p, "--out",  help="Output prefix", default="output/eigs_rarpack")
argv <- parse_args(p)

# Load genotype matrix
cat("Loading genotype matrix from:", argv$geno, "\n")
X <- as.matrix(read.table(argv$geno, header=TRUE))
n <- nrow(X)
m <- ncol(X)
cat("Matrix dimensions:", n, "x", m, "\n")

# Center columns (SNP means)
X1 <- sweep(X, 2, colMeans(X))

# Load A_min
cat("Loading A_min from:", argv$amin, "\n")
A_min <- readRDS(argv$amin)
cat("A_min =", round(A_min, 6), "\n")

# Define d
d <- rowMeans(X1^2 - 1) / A_min

# Load Theta_prod operator
source(file.path("scripts", "Phi_prod.R"))

k <- min(as.integer(argv$k), n - 1)
args_op <- list(X1 = X1, A_min = A_min, d = d)

t0 <- Sys.time()
eigs <- eigs_sym(Theta_prod, k = k, n = n, args = args_op)
runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# Save outputs
dir.create(dirname(argv$out), recursive=TRUE, showWarnings=FALSE)

write.table(eigs$values, paste0(argv$out, "_values.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(eigs$vectors, paste0(argv$out, "_vectors.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# runtime 
write.csv(data.frame(runtime_sec = runtime),
          paste0(argv$out, "_runtime.csv"),
          row.names = FALSE)

cat("rARPACK runtime (sec):", runtime, "\n")





