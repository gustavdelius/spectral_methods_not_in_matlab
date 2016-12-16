
#periodic spectral differentiationn (via matrices)

library("Matrix")


N <- 24
h <- 2*pi/N
x <- h*(1:N)


# I could not find a short way of doing the same type of toeplitz command 
# that matlab uses so I had to code it in a more lengthly way.
# GWD: the following was inspired by p.69 in The Art of R Programming

coluu <- c(0,(.5*(-1)^(1:(N-1)))*1/tan(1:(N-1)*h/2))
m <- matrix(data = 0, nrow = N, ncol = N)
m <- matrix(coluu[(row(m)-col(m)) %% N + 1], nrow=N)

#Differentiation of a hat function

# I imagine there is a faster way to do this `keep possitive values`
# operation. In matlab it can be vectorized

#v <- sapply((1-abs(x-pi)/2), function(v) max(v,0))

# I did Gustav's vectorized version of this command instead

v <- 1 - abs(x - pi) / 2
v[v < 0] <- 0

plot(x,v, type = "b")

plot(x,m %*% v, type = "b", xlab = "x", ylab = "Derivative")

#Differentiation of exp(sin(x))

v <- exp(sin(x))
vprime <- cos(x)*v
plot(x,v, type = "b")
plot(x,m %*% v, type = "b", xlab = "x", ylab = "Derivative")

#compute norm inf error of spectral differentiation

max(abs(m %*% v - vprime))

