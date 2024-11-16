import numpy as np
import matplotlib.pyplot as plt

N = 1000
L = 5
x = np.linspace(-L, L, N+1)
dx = x[1] - x[0]
w = 0.1

def integral(f, dx):
    return np.sum(f*dx, axis=0)

T = -1 / (2 * dx**2) * (np.diag(-2*np.ones(N-1)) + np.diag(np.ones(N-2),1) + np.diag(np.ones(N-2),-1))

V_flat = w * x[1:-1]**2 / 2
V = np.diag(V_flat)

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

""" Evolution of first two eigenstates"""
coeff_0 = np.zeros_like(x[1:-1])
coeff_0[0] = 1 / np.sqrt(2)
coeff_0[1] = 1 / np.sqrt(2)


""" gaussian packet
sigma = 0.1
k0 = 0.01
psi_0 = np.exp(-(x[1:-1])**2 / (2 * sigma**2)) * np.exp(1j * k0 * x[1:-1])
norm  = integral(np.abs(psi_0)**2, dx)
psi_0 = psi_0 / np.sqrt(norm)

coeff_0 = np.zeros_like(eigenstates[0], dtype=complex)
for j in range(0, N-1):
    coeff_0[j] = integral(np.conj(eigenstates[j]) * psi_0, dx)
"""

timespan = np.linspace(0, 100, 1000)

for t in timespan:
    c_n = coeff_0 * np.exp(-1j*En*t)
    psi = eigenstates.T @ (c_n)
    energy = np.sum(En * np.abs(c_n)**2)

    plt.clf()
    #for es in eigenstates[:4]:
    #    plt.plot(x[1:-1], np.abs(es)**2)

    plt.plot(x[1:-1], V_flat)
    plt.plot(x[1:-1], np.abs(psi)**2)
    plt.ylim((0, 2))
    plt.annotate(f'Energy: {energy:.4f}', (-5, 2))
    plt.draw()
    plt.pause(0.01)
