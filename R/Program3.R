# band-limited interpolation

h=1
xmax=10
x = seq(-xmax, xmax, by = h)               # computational grid
xx = seq(-xmax-h/20,xmax+h/20, by = h/10)  # plotting grid

par(mfcol=c(3,1), mar=c(2.5, 4.1, 1, 2.1))
for (plt in 1:3) {
    switch(plt,
           v <- (x==0),                    # delta function
           v <- (abs(x)<=3),               # square wave
           {v <- 1-abs(x)/3; v[v<0]<-0}    # hat function
           )

    p <- rep(0,length(xx));
    for (i in 1:(length(v))) {
        p <- p + v[i]*sin(pi*(xx-x[i])/h)/(pi*(xx-x[i])/h)
    }

    plot(xx, p, type="l", xlab="", ylab="", ylim=c(-.5, 1.5), 
         yaxp  = c(0, 1, 1), las=1)
    points(x, v, col="blue")
    abline(h=1, lty="dotted")
    abline(h=0, lty="dotted")
}
