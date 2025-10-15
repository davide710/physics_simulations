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
while t < 10:
    xs = Xs.copy()
    xs[:, 0] += t #* Xs[:, 1]
    
    plt.clf()
    plt.scatter(xs[:, 0], xs[:, 1])
    plt.xlim(0, L)
    plt.ylim(0, L)
    plt.draw()
    plt.pause(0.05)
    t += 0.1