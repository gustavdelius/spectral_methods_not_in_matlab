# p2.m - convergence of periodic spectral method (compare p1.m)

from scipy.linalg import toeplitz
from numpy import pi, arange, exp, sin, cos, zeros, tan, inf
from numpy.linalg import norm
from matplotlib.pyplot import figure, loglog, hold, grid, xlabel, ylabel, title

figure(figsize=(10, 5))

# Run the following with a range of N
for N in range(2, 100, 2):
    # Set up grid in [-pi,pi] and function u(x)
    h = 2.0*pi/N
    x = -pi + arange(1, N+1)*h
    u = exp(sin(x))
    uprime = cos(x)*u
    # Construct spectral differentiation matrix:
    col = zeros(N)
    col[1:] = 0.5*(-1.0)**arange(1, N)/tan(arange(1, N)*h/2.0)
    row = zeros(N)
    row[0] = col[0]
    row[1:] = col[N-1:0:-1]
    D = toeplitz(col, row)
    # Plot max(abs(D*u-uprime)):
    error = norm(D.dot(u)-uprime, ord=inf)
    loglog(N, error, 'or')
    hold(True)

grid(True, which='both')
xlabel('N')
ylabel('error')
title('Convergence of spectral differentiation')
