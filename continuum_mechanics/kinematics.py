import numpy as np
import matplotlib.pyplot as plt


N = 1e4
L = 10

# position N points in a 2d grid from -L to L, -L to L
Xs = []
for i in range(int(np.sqrt(N))):
    for j in range(int(np.sqrt(N))):
        Xs.append([- L + i * 2 * L / (np.sqrt(N) - 1), -L + j * 2 * L / (np.sqrt(N) - 1)])
Xs = np.array(Xs)

t = 0
xs = Xs.copy()
dt = 0.1
while t < 10:
    xs_spatial = xs.copy()
    xs_spatial[:, 0] += xs[:, 0] * dt / (1 + t)
    xs_spatial[:, 1] += 2 * xs[:, 1] * dt / (1 + t)

    xs[:, 0] = Xs[:, 0] * (1 + t)
    xs[:, 1] = Xs[:, 1] * (1 + t)**2

    
    plt.clf()
    plt.scatter(xs[:, 0], xs[:, 1])
    plt.scatter(xs_spatial[:, 0], xs_spatial[:, 1], color='red', marker='x')
    plt.xlim(-L, L)
    plt.ylim(-L, L)
    plt.draw()
    plt.pause(0.05)
    t += dt