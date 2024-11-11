import numpy as np
import matplotlib.pyplot as plt


timespan = np.linspace(0, 1250, 500)
N = 750
L = 750
x0 = 300
k0 = 0.4
sigma = 10
a = 100
w = 500
V0 = 100

x = np.linspace(0, L, N+1)
dx = x[1] - x[0]

def integral(f):
    return np.sum(f*dx, axis=0)

psi0 = np.exp(-(x[1:-1]-x0)**2 / (2 * sigma**2)) * np.exp(1j * k0 * x[1:-1])
norm  = integral(np.abs(psi0)**2)
psi0 = psi0 / np.sqrt(norm)

#plt.plot(x[1:-1], np.abs(psi0)**2)
#plt.show()

T = -1 / (2 * dx**2) * (np.diag(-2*np.ones(N-1)) + np.diag(np.ones(N-2),1) + np.diag(np.ones(N-2),-1))

V_flat = np.array([V0 if pos < a or pos > a+w else 0 for pos in x[1:-1]])
V = np.diag(V_flat)

En, eigenstates = np.linalg.eigh(T + V)
eigenstates = eigenstates.T
norm = integral(np.abs(eigenstates)**2)
eigenstates = eigenstates / np.sqrt(norm)

c_n0 = np.zeros_like(eigenstates[0], dtype=complex)
for j in range(0, N-1):
    c_n0[j] = integral(np.conj(eigenstates[j]) * psi0)


for t in timespan:
    c_n = c_n0 * np.exp(-1j*En*t)
    psi = eigenstates.T @ (c_n)
    energy = np.sum(En * np.abs(c_n)**2)

    plt.clf()
    for es in eigenstates[:4]:
        plt.plot(x[1:-1], np.abs(es)**2)

    plt.plot(x[1:-1], V_flat)
    plt.plot(x[1:-1], np.abs(psi)**2)
    plt.ylim((0, 0.05))
    plt.annotate(f'Inside: {np.sqrt(integral(np.abs(psi[a:a+w])**2)):.2f}', (0, 0.04))
    plt.annotate(f'Energy: {energy:.4f}', (0, 0.03))
    plt.draw()
    plt.pause(0.01)

