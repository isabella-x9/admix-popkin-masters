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

## References


