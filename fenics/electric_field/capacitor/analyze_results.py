import numpy as np
from scipy.interpolate import LinearNDInterpolator


data = np.loadtxt('electric_field/capacitor/results.txt', delimiter=',')

x, y, z, V = data[:,0], data[:,1], data[:,2], data[:,3]

interp = LinearNDInterpolator(list(zip(x, y, z)), V)

print(f'V(5, 5, 0.1) = {interp(np.array([5, 5, 0.15]))}')
print(f'V(5, 5, 0.3) = {interp(np.array([5, 5, 0.25]))}')
