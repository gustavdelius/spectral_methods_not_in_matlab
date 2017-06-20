# Convergence of periodic spectral method (compare to Program1.R)

DerivativeMatrixErrorInf <- function(N) {
    # set up grid in [-pi,pi] and function u(x):
    h <- 2*pi/N
    x <- -pi+(1:N)*h
    u <- exp(sin(x))
    uprime <- cos(x)*u
    
    # Construct spectral differentiation matrix:
    column <- c(0, .5*(-1)^(1:(N-1))/tan(1:(N-1)*h/2))
    m <- matrix(data = 0, nrow = N, ncol = N)
    for (i in 1:N){
        for (j in 1:N){
            m[i,j] <- column[((i-j) %% N)+1]
        }   
    }
    return(max(abs((m %*% u) - uprime)))  
}


Nvec <- 2*(1:50)
errorvalues <- sapply(2*Nvec, DerivativeMatrixErrorInf)

plot(Nvec, errorvalues, log="xy", xlab="N", ylab="error",
     main="Convergence of spectral differentiation")
