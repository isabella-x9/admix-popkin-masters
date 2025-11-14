# Define function Phi * v 

Phi_prod <- function(x, args) {
  X1 <- args$X1
  A_min <- args$A_min
  m <- ncol(X1)
  
  term1 <- sum(x) * (A_min + 1 )
  term2 <- (1 / m) * (X1 %*% crossprod(X1, x))
  ( term1 - term2 ) / A_min 
}
