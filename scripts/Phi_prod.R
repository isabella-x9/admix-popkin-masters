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
