# Define function Phi * v 

Phi_prod <- function(x, args) {
  X1 <- args$X1
  m <- ncol(X1)
  A_min <- args$A_min
  
  term1 <- sum(x) * (A_min + 1)
  term2 <- (1 / m) * (X1 %*% crossprod(X1, x))
  ( term1 - term2 ) / A_min 
}


# Define function Theta product 
Theta_prod <- function(x, args) {
  phi_x <- Phi_prod(x, args)
  phi_x - args$d * x
}


### Define a new version for handling missingness
Phi_prod_miss <- function(x, args) {
  p_obs <- args$p_obs
  X1 <- args$X1
  m <- ncol(X1)
  A_min <- args$A_min
  
  term1 <- sum(x) * (A_min + 1)
  term2 <- (1 / m) * drop(X1 %*% crossprod(X1, x / p_obs)) / p_obs
  ( term1 - term2 ) / A_min 
}

### Define a new version for handling missingness
Theta_prod_miss <- function(x, args) {
  phi_x <- Phi_prod_miss(x, args)
  phi_x - args$d * x
}
