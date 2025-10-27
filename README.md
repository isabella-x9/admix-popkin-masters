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
├── scripts/
│   ├── 01_build_popkin_grm.R
│   ├── 02_popkin_evd.R
│
├── data/
│   ├── test_geno.tsv
│   └── .gitkeep
│
├── output/
│   └── .gitkeep
│
└── README.md
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
- Calculates marker weights based on allele frequency variance
- Saves outputs to:
  - `output/Phi.tsv`
  - `output/Phi_marker_weights.tsv`

### Stage 02 - Popkin Eigendecomposition (EVD)
- Loads $\Phi$ from Stage 01 and performs a **top-K eigendecomposition using** `RSpectra::eigs_sym()`
- Produces leading eigenvalues and eigenvectors for downstream analysis
- Saves results to: 
  - `output/eigen_vals.tsv`
  - `output/eigen_vecs.tsv`

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

Expected outputs appear in the `outputs/` folder. 

---

## References


