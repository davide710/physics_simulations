import numpy as np
import matplotlib.pyplot as plt

N = 1000
L = 5
x = np.linspace(-L, L, N+1)
dx = x[1] - x[0]
w = 1

def integral(f, dx):
    return np.sum(f*dx, axis=0)

T = -1 / (2 * dx**2) * (np.diag(-2*np.ones(N-1)) + np.diag(np.ones(N-2),1) + np.diag(np.ones(N-2),-1))

V = w * x[1:-1]**2 / 2
V = np.diag(V)

H = T + V
En, eigenstates = np.linalg.eigh(T + V)
eigenstates = eigenstates.T
norm = integral(np.abs(eigenstates)**2, dx)
eigenstates = eigenstates / np.sqrt(norm)

""" the first eigenstate is 1/pi^(1/4) exp(-x^2/2) as it should
plt.plot(x[1:-1], np.abs(eigenstates[0])**2, label=f'computed, {integral(np.abs(eigenstates[0])**2, dx)}')
plt.plot(x, np.abs(1 / np.power(np.pi, 1/4) * np.exp(-x**2 / 2))**2, label=f'exact, {integral(np.abs(1 / np.power(np.pi, 1/4) * np.exp(-x**2 / 2))**2, dx)}')
plt.legend()
plt.show()
"""

