# p5.m - repetition of p4.m via FFT

# For complex v, delete "real" commands.

from numpy import pi, maximum, abs, zeros, arange, real, sin, cos, exp, inf
from numpy.fft import fft, ifft
from numpy.linalg import norm
from matplotlib.pyplot import figure, subplot, plot, axis, title, text

figure(figsize=(10, 6))

# Set up grid:
N = 24
h = 2*pi/N
x = h*arange(1, N+1)

# Differentiation of a hat function:
v = maximum(0, 1-abs(x-pi)/2)
v_hat = fft(v)
w_hat = 1j*zeros(N)
w_hat[0:N//2] = 1j*arange(0, N//2)
w_hat[N//2+1:] = 1j*arange(-N//2+1, 0, 1)
w_hat = w_hat * v_hat
w = real(ifft(w_hat))

subplot(3, 2, 1)
plot(x, v, '.-')
axis([0, 2*pi, -.5, 1.5])
title('function')
subplot(3, 2, 2)
plot(x, w, '.-')
axis([0, 2*pi, -1, 1])
title('spectral derivative')

# Differentiation of exp(sin(x)):
v = exp(sin(x))
vprime = cos(x)*v
v_hat = fft(v)
w_hat = 1j*zeros(N)
w_hat[0:N//2] = 1j*arange(0, N//2)
w_hat[N//2+1:] = 1j*arange(-N//2+1, 0, 1)
w_hat = w_hat * v_hat
w = real(ifft(w_hat))

subplot(3, 2, 3)
plot(x, v, '.-')
axis([0, 2*pi, 0, 3])
subplot(3, 2, 4)
plot(x, w, '.-')
axis([0, 2*pi, -2, 2])
error = norm(w-vprime, inf)
text(2.0, 1.4, "max error="+str(error))
