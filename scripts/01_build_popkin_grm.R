# 01_build_popkin_grm.R
# ------------------------------------------------
# Description:
#   This script calculates the Popkin kinship matrix (Phi) and marker weights (M)
#   from a general genotype matrix file. 
# ============================================================

suppressPackageStartupMessages({
  library(popkin)
  library(genio)
  library(argparser)
})

# Parse command-line arguments 
p <- arg_parser("Compute Popkin kinship and marker weights from a genotype matrix")
p <- add_argument(p, "--geno", help = "Path to genotype file (tab- or space-delimited, numeric matrix)", default = NULL)
p <- add_argument(p, "--out",  help = "Output prefix for GRM (no extension)", default = NULL)
argv <- parse_args(p)

# manual checks
if (is.null(argv$geno) || is.null(argv$out)) {
  cat("\nERROR: --geno and --out are required.\n\n")
  print(p)
  quit(status = 2)
}

geno_path <- argv$geno
base_out  <- argv$out

# Load genotype data
message("Reading genotype matrix from: ", geno_path)
X <- as.matrix(read.table(geno_path, header = FALSE, sep = "", check.names = FALSE))
storage.mode(X) <- "numeric"
message("Dimensions: ", paste(dim(X), collapse = " × "))

# Compute kinship 
message("Computing Popkin kinship matrix ...")
obj <- popkin(X, want_M = TRUE)
kinship <- obj$kinship
M <- obj$M

# Save results
message("Writing GRM to: ", base_out)
write_grm(base_out, kinship, M = M)
message("RM successfully written.")






