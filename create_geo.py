import numpy as np

N_POINTS = 100

with open('spira.geo', 'w') as f:
    for i in range(N_POINTS):
        f.write('Point('+ str(i+1) + ') = {' + str(10*np.cos(2*np.pi/N_POINTS*i)) + ', ' + str(10*np.sin(2*np.pi/N_POINTS*i)) + ', 0, 1.0};\n//+\n')

with open('spira.geo', 'a') as f:
    for i in range(N_POINTS-1):
        f.write('Line('+ str(i+1) + ') = {' + str(i+1) + ', ' + str(i+2) +'};\n//+\n')
    f.write('Line(' + str(N_POINTS) + ') = {' + str(N_POINTS) + ', 1};\n//+\n')

with open('spira.geo', 'a') as f:
    f.write('Physical Curve("spira") = {' + str([i for i in range(1, N_POINTS+1)]) + '};\n//+\n')

