# Variable coefficient wave equation pde evolution by spectral methods

library("rgl")

# grid setup
N <- 128
h <- 2*pi/N
x <- h*(1:N)
t <- 0
dt <- h/4  # time step of one iteration in de solver
c <- 0.2 + sin(x-1)^2 # variable coefficient
# Initial function
v <- exp(-20*(x-1)^2)
# Function just before starting point (bit of a hack, better to work it out
# using backwards time euler)
vold <- exp(-20*(x-0.2*dt-1)^2)

# leap frog time stepping
tmax <- 8 # total amount of t system is run for
tplot <- 0.15  # approximate amount of time between plotted points
# number of iterations of the solver between plotted points
plotgap <- round(tplot/dt)
# number of time points data is plotted for
nplots <- round(tmax/tplot)

# Create array to store PDE evolution data
# A row stores the data from a given time step,
# the first value on that row is the time, the rest
# are the amplitute values at the N grid points
# this yields a similar output to that given by the ODE function
spaceTimeData <- matrix(data = 0,
                        nrow = nplots + 1,
                        ncol = N + 1)

spaceTimeData[1, ] <- c(0, v)  # First row holds initial time and intial value

#loop for ith plot output
for (i in 1:nplots) {
    #loop for nth time step between plotting moments
    for (n in 1:plotgap) {
        t <- t + dt
        vHat <- fft(v)
        wHat <- 1i * (c(0:(N/2-1), 0, (-N/2+1):-1))*vHat
        w <- Re(fft(wHat, inverse=TRUE)/N)
        # Use leapfrog method to calculate v at next timestep
        vnew <- vold - 2*dt*c*w
        vold <- v
        v <- vnew
    }
    spaceTimeData[i+1, ] <- c(t, v)
}

# I think this code evolves the PDE just like in the book, although it looks
# dodgy because of the leapfrog initialization hack Trefethen used.

# Plot solution at final time
plot(x, spaceTimeData[nplots+1, 2:(N+1)], type="l")

# Make 3D plot of solution over time
Time <- spaceTimeData[, 1]
Amplitude <- spaceTimeData[, 2:(N+1)]
persp3d(Time, x, Amplitude, col = "lightblue")
# I find it strange that the function is plotted with x increasing from right 
# to left
