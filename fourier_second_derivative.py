import numpy as np
import matplotlib.pyplot as plt

N = 1000
L = 10

x = np.linspace(-L/2, L/2, N)
f = 1 / (1 + x**2)

f_xx_exact = 8 * x**2 / (1 + x**2)**3 - 2 / (1 + x**2)**2

f_hat = np.fft.fft(f)
k = 2*np.pi/L * np.arange(-N/2, N/2)
k = np.fft.fftshift(k)
f_xx = np.real(np.fft.ifft(-k**2 * f_hat))

dx = x[1] - x[0]
D2_matrix = 1 / dx**2 * (np.diag(-2*np.ones(N-1)) + np.diag(np.ones(N-2),1) + np.diag(np.ones(N-2),-1))
f_xx_finite_differences = D2_matrix @ f[1:]
f_xx_finite_differences[0] = 0
f_xx_finite_differences[-1] = 0

plt.plot(x, f)
plt.plot(x, f_xx_exact, label='exact')
plt.scatter(x, f_xx, label='fft')
plt.scatter(x[1:], f_xx_finite_differences, label='finite diff')
plt.legend()
plt.show()
