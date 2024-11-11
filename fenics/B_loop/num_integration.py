import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def f(x):
    R = 6.5
    dt = 2*np.pi / 1000
    I = 0
    for i in range(1000):
        t = i*dt
        I += dt * (1 - x*np.cos(t)/R) * (1 + R**2/x**2 - 2*R*np.cos(t)/x)**(-3/2)
    return 4e-7*R*R*I / x**3

fig, axs = plt.subplots(1, 2, constrained_layout=True, sharey='row')

x = np.linspace(0.1, 5, num=1000)
y = np.array([f(u) for u in x])
axs[0].plot(x, y)

df = pd.read_csv('results.csv')
axs[1].plot(df.r, df.B)

plt.show()