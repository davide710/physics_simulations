import numpy as np
import matplotlib.pyplot as plt


timespan = np.linspace(0, 1500, 600)
N = 750
L = 750
x0 = 350
k0 = 0.1
sigma = 10

x = np.linspace(0, L, N+1)
dx = x[1] - x[0]

def integral(f):
    return np.sum(f*dx, axis=0)

psi0 = np.exp(-1/2 * (x[1:-1]-x0)**2 / sigma**2) * np.exp(1j * k0 * x[1:-1])
norm  = integral(np.abs(psi0)**2)
psi0 = psi0 / np.sqrt(norm)

T = -1/2 * 1/dx**2 * (np.diag(-2*np.ones(N-1))+ np.diag(np.ones(N-2),1)+ np.diag(np.ones(N-2),-1))

En, eigenstates = np.linalg.eigh(T)
eigenstates = eigenstates.T
norm = integral(np.abs(eigenstates)**2)
eigenstates = eigenstates / np.sqrt(norm)

c_n = np.zeros_like(eigenstates[0], dtype=complex)
for j in range(0, N-1):
    c_n[j] = integral(np.conj(eigenstates[j]) * psi0)


for t in timespan:
    psi = eigenstates.T @ (c_n * np.exp(-1j*En*t))
    
    plt.clf()
    plt.plot(x[1:-1], np.abs(psi)**2)
    plt.draw()
    plt.pause(0.01)

