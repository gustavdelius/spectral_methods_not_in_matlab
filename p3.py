# p3.m - band-limited interpolation

from numpy import arange, maximum, abs, zeros, sin, pi
from matplotlib.pyplot import subplot, figure, plot, grid, axis

h = 1.0
xmax = 10.0
x = arange(-xmax, xmax+h, h)              # computational grid
xx = arange(-xmax-h/20, xmax+h/20, h/10)  # plotting grid
figure(figsize=(10, 10))
for pl in range(3):
    subplot(4, 1, pl+1)
    if pl == 0:
        v = (x == 0)                      # Kroneker delta
    elif pl == 1:
        v = (abs(x) <= 3.0)               # square wave
    else:
        v = maximum(0.0, 1.0-abs(x)/3.0)  # hat function
    plot(x, v, '.')
    grid(True)
    p = zeros(len(xx))
    for i in range(len(x)):
        p = p + v[i]*sin(pi*(xx-x[i])/h)/(pi*(xx-x[i])/h)
    plot(xx, p)
    axis([-xmax, xmax, -0.5, 1.5])
