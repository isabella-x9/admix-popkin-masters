# Master's Project: Scalable Calculation of the Top K Eigenvectors and Eigenvalues of the Popkin Kinship Matrix Estimate

**Author:** Isabella Xu  
**Advisor:** Dr. Alex Ochoa, Duke University  
**Repository:** [isabella-x9/admix-popkin-masters](https://github.com/isabella-x9/admix-popkin-masters)  

---

## Overview
This project extends computational methods for large-scale population genetics by developing a **scalable eigen-decomposition pipeline** for the **Popkin kinship matrix**, used in **AdmixCor** for admixture inference.  

The goal is to compute the top *K* eigenvectors and eigenvalues efficiently for large genomic datasets, hence enabling downstream analysis of population structure, kinship, and admixture proportions. 

---

## Repository Structure (tentative)
```
admix-popkin-masters/
├── admix-popkin-masters.Rproj
├── README.md
│
├── data/
│   ├── test_geno.tsv
│   └── raw/
|       └── genotypes.tsv 
│
├── scripts/
│   ├── 01_build_popkin_grm.R         # Stage 1: Popkin Phi, marker weights, A_min
│   ├── 02_popkin_evd.R               # Stage 2: RSpectra eigendecomposition
│   ├── 03_eigs_rarpack.R             # Stage 3: rARPACK eigendecomposition
│   ├── 04_runtime_scaling.R          # Stage 4: runtime benchmark
│   └── Phi_prod.R                    # Phi*v matrix–vector product for Lanczos
│
├── output/
│   ├── Phi.tsv                       # Popkin kinship matrix (Stage 1)
│   ├── Phi_marker_weights.tsv        # Marker weights (Stage 1)
│   ├── Phi_Amin.rds                  # A_min scalar used in Popkin-EVD (Stage 1)
│   │
│   ├── eigen_vals.tsv                # Eigenvalues from explicit Phi (Stage 2)
│   ├── eigen_vecs.tsv                # Eigenvectors from explicit Phi (Stage 2)
│   │
│   ├── rarpack_vals.tsv              # Eigenvalues from rARPACK (Stage 3)
│   ├── rarpack_vecs.tsv              # Eigenvectors from rARPACK (Stage 3)
│   ├── rarpack_runtime.rds           # Runtime for rARPACK EVD (Stage 3)
│   │
│   ├── runtime_scaling.csv           # RSpectra vs. rARPACK runtime results (Stage 4)
│   ├── runtime_scaling.png           # Plot of runtime scaling (Stage 4)
│   │
│   └── grm/                          # Additional GRM files (optional/experimental)
│       ├── test_popkin.grm.bin
│       ├── test_popkin.grm.id
│       └── test_popkin.grm.N.bin
│
├── renv/
├── renv.lock
└── .gitignore
```

---

## Setup and Reproducibility
To reproduce an environment used for this project: 

```r
# Clone the repository
git clone git@github.com:isabella-x9/admix-popkin-masters.git
cd admix-popkin-masters

# Restore R packages
install.packages("renv")
renv::restore()
```

---

## Current Workflow
### Stage 01 – Build Popkin GRM
- Loads a small genotype matrix (`data/test_geno.tsv`)
- Computes the **Popkin kinship matrix** using the `popkin` R package
- Calculates marker weights and the required `A_min`
- Saves outputs to:
  - `output/Phi.tsv`
  - `output/Phi_marker_weights.tsv`
  - `output/Phi_Amin.rds`

### Stage 02 - Popkin Eigendecomposition (EVD)
- Loads $\hat{\Phi}$ from Stage 01 and P
- Performs a **top-K eigendecomposition using** `RSpectra::eigs_sym()`
- Produces leading eigenvalues and eigenvectors for downstream analysis
- Saves results to: 
  - `output/eigen_vals.tsv`
  - `output/eigen_vecs.tsv`

### Stage 03 - Scalable Eigendecomposition (`rARPACK`)
- Avoids forming $\hat{\Phi}$ explicitly
- Uses a custom `Phi_prod.R` function that computes $\Phi \cdot v$ via the Popkin-EVD formulation
- Uses `rARPACK::eigs_sym()` for scalable eigendecomposition
- Produces identical eigenvalues to RSpectra after correcting `X1 = X - 1`
- Saves results to: 
  - `output/rarpack_vals.tsv`
  - `output/rarpack_vecs.tsv`
  - `output/rarpack_runtime.rds`

### Stage 04 – Runtime Scaling
- Compares runtime for: 
  - `RSpectra` (explicit matrix)
  - `rARPACK` (matrix-free)
- Saves: 
  - `output/runtime_scaling.csv`
  - `output/runtime_scaling.png` 
---

## Quick Test Run

To verify the pipeline end-to-end with toy data:

```r
# Create small random genotype dataset
dir.create("data", showWarnings = FALSE)
set.seed(1)
X <- matrix(sample(0:2, 10*5, replace = TRUE, prob = c(0.49, 0.42, 0.09)),
            nrow = 10, ncol = 5)
write.table(X, "data/test_geno.tsv", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)

# Stage 01: Build kinship matrix
source("scripts/01_build_popkin_grm.R")

# Stage 02: Eigendecomposition
source("scripts/02_popkin_evd.R")
```

Expected outputs appear in the `output/` folder. 

---

## References

- Ochoa & Storey (2021). *Estimating FST and kinship for arbitrary population structures*. PLoS Genetics. 
- Ochoa et al. (2025). *AdmixCor: Admixture inference from genetic covariance*. Duke University, Ochoa Lab.
- Popkin R vignette: https://cran.r-project.org/web/packages/popkin/vignettes/popkin.html

