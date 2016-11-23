library("Matrix")




DerivativeMatrixError <- function(N){
  h <- 2*pi/N;
  x <- -pi+(1:N)*h;
  u <- exp(sin(x));
  uprime <- cos(x)*u;
  D <- sparseMatrix((1:N), c((2:N),1), x = 2*rep(1, N)/3)-sparseMatrix((1:N), c((3:N),1,2), x = rep(1, N)/12)
  D <- (D-t(D))/h;
return(max(abs(D %*% u -uprime))
)  
}


Nvec=2^(3:12);
errorvalues <- rep(0,length(Nvec));
for (i in 1:(length(errorvalues))) {
 errorvalues[i] <- DerivativeMatrixError(Nvec[i])
}

plot(Nvec,errorvalues, log="xy")
