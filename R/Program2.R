library("Matrix")

DerivativeMatrixErrorInf <- function(N){
    h <- 2*pi/N;
    x <- -pi+(1:N)*h;
    u <- exp(sin(x));
    uprime <- cos(x)*u;
    
    coluu=c(0,(.5*(-1)^(1:(N-1)))*1/tan(1:(N-1)*h/2))
    m=matrix(data = 0, nrow = N, ncol = N)
    for (i in 1:N){
        for (j in 1:N){
            m[i,j]=coluu[((i-j) %% N)+1]
        }   
    }
    return(max(abs((m %*% u) -uprime)))  
}

Nvec=(1:50);
errorvalues2 <- sapply(2*Nvec, DerivativeMatrixErrorInf)

plot(2*Nvec,errorvalues2, log="xy")