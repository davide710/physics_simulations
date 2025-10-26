import numpy as np
import matplotlib.pyplot as plt


N = 1e2
L = 10

# position N points in a 2d grid from -L to L, -L to L
Xs = []
for i in range(int(np.sqrt(N))):
    for j in range(int(np.sqrt(N))):
        Xs.append([- L + i * 2 * L / (np.sqrt(N) - 1), -L + j * 2 * L / (np.sqrt(N) - 1)])
Xs = np.array(Xs)

t = 0
xs = Xs.copy()
xs_spatial = Xs.copy()
dt = 0.1
while t < 10:
    xs_spatial[:, 0] += -2 * xs_spatial[:, 1] * dt
    xs_spatial[:, 1] += 2 * xs_spatial[:, 0] * dt

    xs[:, 0] = -Xs[:, 1] * np.sin(2*t) + Xs[:, 0] * np.cos(2*t)
    xs[:, 1] = Xs[:, 1] * np.cos(2*t) + Xs[:, 0] * np.sin(2*t)

    
    plt.clf()
    plt.scatter(xs[:, 0], xs[:, 1])
    plt.scatter(xs_spatial[:, 0], xs_spatial[:, 1], color='red', marker='x')
    plt.xlim(-L, L)
    plt.ylim(-L, L)
    plt.draw()
    plt.pause(0.05)
    t += dt