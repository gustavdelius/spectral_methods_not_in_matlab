

#variable coefficient wave equation pde
# evolution by spectral methods

library("rgl")

# grid setup

N <- 128
h <- 2*pi/N
x <- h*(1:N)
t <- 0

#dt = time step of one iteration in de solver

dt <- h/4

#variable coefficient

c <- 0.2 + sin(x-1)^2

#Initial funtction

v <- exp(-100*(x-1)^2)

#Function just before starting point (bit of a hack, better to work it out using
#backwards time euler)

vold <- exp(-100*(x-0.2*dt-1)^2)

#leap frog time stepping

#tmax is total amount of t system is run for, as on plot

tmax <- 8

#tplot = ?amount of time between plotted points ?

tplot <- 0.15

#plotgap = number of iterations of de solver between plotted points

plotgap <- round(tplot/dt)

#nplots = number of time points data is plotted for

nplots <- round(tmax/tplot)


#Use row binding to store PDE evolution data
# A row stores the data from a given time step, 
# the first value on that row is the time, the rest 
# are the amplitute values at the N grid points
# this yields a similar output to that given by the ODE function 
#we use rbind() to keep concatenating new rows of generated values

spaceTimeData <- matrix(data = 0, nrow = nplots + 1, ncol = N+1)


spaceTimeData <- c(0,v)

#loop for ith plot output
for (i in 1:nplots){
    #loop for nth time step between plotting moments
        for (n in 1:plotgap){
            t <- t+dt
            
            vHat <- fft(v)
            wHat <- 1i*(c(0:(N/2-1),0,(-N/2+1):-1))*vHat
            w <- Re(fft(wHat, inverse=TRUE)/N)
            
            vnew <- vold - 2*dt*c*w
            vold <- v
            v <- vnew
        } 
    spaceTimeData <- rbind(spaceTimeData, c(t,v) )
}

# I think this code evolves the PDE just like in the book, although it looks 
# dodgy because of the leapfrog initialization hack Trefethen used.

plot(x,spaceTimeData[nplots+1,2:(N+1)])

Time <- spaceTimeData[,1]

Amplitude <- spaceTimeData[,(2:(N+1))]

# I find it strange that the function is plotted with x increasing from right to left


persp3d(Time, x, Amplitude, col = "lightblue")


