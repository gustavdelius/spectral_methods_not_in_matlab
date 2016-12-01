library("Matrix")


h=1
xmax=10
x = seq(-xmax, xmax, by = h)
xx = seq(-xmax-h/20,xmax+h/20, by = h/10)

#Interpolation using the delta function
v = (x==0)*1

p <- rep(0,length(xx));
for (i in 1:(length(v))) {
    p = p + v[i]*sin(pi*(xx-x[i])/h)/(pi*(xx-x[i]))
}

plot(x,v,col="blue")
lines(xx,p,col="green")

#Interpolation using the square wave function
v = (abs(x)<=3)*1

p <- rep(0,length(xx));
for (i in 1:(length(v))) {
    p = p + v[i]*sin(pi*(xx-x[i])/h)/(pi*(xx-x[i]))
}

plot(x,v,col="blue")
lines(xx,p,col="green")


#Interpolation using the hat function
v = 1-abs(x)/3
v = v*(v>=0)

p <- rep(0,length(xx));
for (i in 1:(length(v))) {
    p = p + v[i]*sin(pi*(xx-x[i])/h)/(pi*(xx-x[i]))
}

plot(x,v,col="blue")
lines(xx,p,col="green")

