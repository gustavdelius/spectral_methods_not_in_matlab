# p1.m - convergence of fourth-order finite differences

from scipy.sparse import coo_matrix
from numpy import arange, pi, exp, sin, cos, ones, inf
from numpy.linalg import norm
from matplotlib.pyplot import figure, loglog, hold, semilogy, text, grid, \
                              xlabel, ylabel, title

figure(figsize=(10, 5))

# Run the following with a range of N
Nvec = 2**arange(3, 13)
for N in Nvec:
    # Set up grid in [-pi,pi] and function u(x)
    h = 2*pi/N
    x = -pi + arange(1, N+1)*h
    u = exp(sin(x))
    uprime = cos(x)*u
    # Construct sparse 4th-order differentiation matrix:
    e = ones(N)
    e1 = arange(0, N)
    e2 = arange(1, N+1)
    e2[N-1] = 0
    e3 = arange(2, N+2)
    e3[N-2] = 0
    e3[N-1] = 1
    D = coo_matrix((2*e/3, (e1, e2)), shape=(N, N)) \
        - coo_matrix((e/12, (e1, e3)), shape=(N, N))
    D = (D - D.T)/h
    error = norm(D.dot(u)-uprime, ord=inf)
    loglog(N, error, 'or')
    hold(True)

# Add line N^-4 to plot
semilogy(Nvec, Nvec**(-4.0), '--')
text(105, 5e-8, '$N^{-4}$', fontsize=20)
grid(True, which='both')
xlabel('N')
ylabel('error')
title('Convergence of fourth-order finite difference')
