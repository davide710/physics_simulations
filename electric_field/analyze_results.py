import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv('results_x.csv')

rho = 1e-3
eps = 8.85e-12
r_int = 0.05
h = 1
q = 3.14 * r_int**2 * h * rho

r = df.r

V_sim = df.V

V_pred = q * np.log(r/0.05) / 2 / eps / 3.14 / h

plt.scatter(r, -V_sim)
plt.plot(r, V_pred)
plt.show()

