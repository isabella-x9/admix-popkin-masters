# 02_popkin_evd.R
# ------------------------------------------------
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
