# convergence of 4th order finite difference method

library("Matrix")

DerivativeMatrixError <- function(N) {
  # set up grid in [-pi,pi] and function u(x):
  h <- 2*pi/N
  x <- -pi+(1:N)*h
  u <- exp(sin(x))
  uprime <- cos(x)*u
  # Construct sparse fourth-order differentiation matrix:
  e <- rep(1, N)
  D <- sparseMatrix(1:N, c(2:N, 1), x = 2*e/3) -
      sparseMatrix(1:N, c(3:N, 1, 2), x = e/12)
  D <- (D-t(D))/h
  return(max(abs(D %*% u - uprime)))  
}


Nvec=2^(3:12)
errorvalues <- sapply(Nvec, DerivativeMatrixError)

plot(Nvec, errorvalues, log="xy", xlab="N", ylab="error",
     main="Convergence of fourth-order finite differences")
lines(Nvec, Nvec^(-4), lty="dashed")
text(105, 5e-8, expression(N^{-4}))
