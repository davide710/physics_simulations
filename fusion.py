import numpy as np
import matplotlib.pyplot as plt

# https://www.astro.utoronto.ca/~mahajan/notebooks/quantum_tunnelling.html

N = 2500
L = 250
x = np.linspace(20, L, N+1)
dx = x[1] - x[0]
a = 40
V0 = 40

def integral(f, dx):
    return np.sum(f*dx, axis=0)

def gaussian_packet(sigma=5, x0=90, k0=-1):
    psi_0 = np.exp(-(x[1:-1] - x0)**2 / (2 * sigma**2)) * np.exp(1j * k0 * x[1:-1])
    norm  = integral(np.abs(psi_0)**2, dx)
    return psi_0 / np.sqrt(norm)

def hamiltonian(V0=40):
    T = -1 / (2 * dx**2) * (np.diag(-2*np.ones(N-1)) + np.diag(np.ones(N-2),1) + np.diag(np.ones(N-2),-1))
    V_flat = np.array([0 if pos < a else V0/pos for pos in x[1:-1]])
    V = np.diag(V_flat)
    H = T + V
    return H, V_flat

def get_eigenstates(H):
    En, eigenstates = np.linalg.eigh(H)
    eigenstates = eigenstates.T
    norm = integral(np.abs(eigenstates)**2, dx)
    eigenstates_list = eigenstates / np.sqrt(norm)
    eigenstates_matrix = eigenstates_list.T
    return En, eigenstates_list, eigenstates_matrix

def get_coeffs_in_basis(psi, basis_list):
    coeffs = np.zeros_like(basis_list[0], dtype=complex)
    for j in range(0, len(basis_list)):
        coeffs[j] = integral(np.conj(basis_list[j]) * psi, dx)
    return coeffs


H, V_flat = hamiltonian(V0=V0)

psi_0 = gaussian_packet(sigma=5, x0=90, k0=-1)

En, eigenstates_list, eigenstates_matrix = get_eigenstates(H)

coeff_0 = get_coeffs_in_basis(psi_0, eigenstates_list)

timespan = np.linspace(0, 100, 100)
for t in timespan:
    c_n = coeff_0 * np.exp(-1j*En*t)
    psi = eigenstates_matrix @ (c_n)
    energy = np.sum(En * np.abs(c_n)**2)
    tunnelling_prob = integral(np.abs(psi[x[1:-1] < a])**2, dx)

    plt.clf()
    plt.plot(x[1:-1], V_flat)
    plt.plot(x[1:-1], np.abs(psi)**2)
    plt.annotate(f'Energy = {np.sum(En * np.abs(c_n)**2):.2f}', (a, 0))
    plt.annotate(f'Tunnelling prob: {tunnelling_prob:.2f}', (a, V0 / a * 0.9))
    plt.draw()
    if tunnelling_prob > 0.2:
        plt.pause(5)
        break
    plt.pause(0.01)
