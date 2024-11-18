import numpy as np
import matplotlib.pyplot as plt

N = 1000
L = 500
x = np.linspace(-L, L, N)
dx = x[1] - x[0]

def integral(f, dx):
    return np.sum(f*dx, axis=0)

T = -1 / (2 * dx**2) * (np.diag(-2*np.ones(N)) + np.diag(np.ones(N-1),1) + np.diag(np.ones(N-1),-1))

H = T
En, eigenstates = np.linalg.eigh(T)
eigenstates = eigenstates.T
norm = integral(np.abs(eigenstates)**2, dx)
eigenstates = eigenstates / np.sqrt(norm)

#plt.plot(x, np.abs(eigenstates[0])**2, label=f'first')
#plt.plot(x, np.abs(eigenstates[1])**2, label=f'second')
#plt.legend()
#plt.show()

sigma = 10
k0 = 100
psi_0 = np.exp(-(x)**2 / (2 * sigma**2)) * np.exp(1j * k0 * x)
norm  = integral(np.abs(psi_0)**2, dx)
psi_0 = psi_0 / np.sqrt(norm)

coeff_0 = np.zeros_like(eigenstates[0], dtype=complex)
for j in range(0, N-1):
    coeff_0[j] = integral(np.conj(eigenstates[j]) * psi_0, dx)

timespan = np.linspace(0, 10, 1000)

for t in timespan:
    c_n = coeff_0 * np.exp(-1j*En*t)
    psi = eigenstates.T @ (c_n)
    energy = np.sum(En * np.abs(c_n)**2)

    plt.clf()

    plt.plot(x, np.abs(psi)**2)
    plt.annotate(f'{x[np.argmax(np.abs(psi)**2)]:.2f}', (0, 0))
    plt.annotate(f'Energy = {np.sum(En * np.abs(c_n)**2):.2f}', (0, 0.3))
    plt.draw()
    plt.pause(0.01)