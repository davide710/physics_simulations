import numpy as np
import matplotlib.pyplot as plt

def f(x):
    dt = 2*np.pi / 1000
    I = 0
    for i in range(1000):
        t = i*dt
        I += dt * (1 - x*np.cos(t)/6) * (1 + 36/x**2 - 12*np.cos(t)/x)**(-3/2)
    return 4e-7*36*I / x**3

x = np.linspace(0.1, 5, num=1000)
y = np.array([f(u) for u in x])
plt.plot(x, y)
plt.show()