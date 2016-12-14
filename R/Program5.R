
#periodic spectral differentiationn (via fft)

library("Matrix")


N <- 24
h <- 2*pi/N
x <- h*(1:N)

#Differentiation of a hat function


#v <- sapply((1-abs(x-pi)/2), function(v) max(v,0))

# I did Gustav's vectorized version of this command instead

v <- 1 - abs(x - pi) / 2
v[v < 0] <- 0

vHat <- fft(v)
wHat <- 1i*(c(0:(N/2-1),0,(-N/2+1):-1))*vHat
w <- Re(fft(wHat, inverse=TRUE)/N)

plot(x,v, type = "b")


plot(x,w, type = "b")


#Differentiation of exp(sin(x))

v <- exp(sin(x))
vprime <- cos(x)*v

vHat <- fft(v)
wHat <- 1i*(c(0:(N/2-1),0,(-N/2+1):-1))*vHat
w <- Re(fft(wHat, inverse=TRUE)/N)

plot(x,v, type = "b")


plot(x,w, type = "b")


#compute norm inf error of spectral differentiation

max(abs(w - vprime))

# A tiny bit more accurate than what matlab does
