import numpy as np

######################### INNER CYLINDER


N_POINTS = 100
R = 0.05
"""
with open('wire.geo', 'w') as f:
    for i in range(N_POINTS):
        f.write('Point('+ str(i+1) + ') = {' + f'{R*np.cos(2*np.pi/N_POINTS*i):.4f}' + ', ' + f'{R*np.sin(2*np.pi/N_POINTS*i):.4f}' + ', 0, 1.0};\n//+\n')
"""
with open('wire.geo', 'a') as f:
    for i in range(N_POINTS, 2 * N_POINTS):
        f.write('Point('+ str(i+1) + ') = {' + f'{R*np.cos(2*np.pi/N_POINTS*i):.4f}' + ', ' + f'{R*np.sin(2*np.pi/N_POINTS*i):.4f}' + ', 10.0, 1.0};\n//+\n')
"""
with open('wire.geo', 'a') as f:
    for i in range(2*N_POINTS-1):
        f.write('Line('+ str(i+1) + ') = {' + str(i+1) + ', ' + str(i+2) +'};\n//+\n')  ####### poi da modificare
    f.write('Line(' + str(N_POINTS) + ') = {' + str(N_POINTS) + ', 1};\n//+\n') ####### poi da modificare

with open('wire.geo', 'a') as f:
    for i in range(N_POINTS):
        f.write('Line('+ str(200 + i + 1) + ') = {' + str(i+1) + ', ' + str(100 + i + 1) +'};\n//+\n')

with open('wire.geo', 'a') as f:
    for i in range(N_POINTS):
        f.write('Curve Loop('+ str(2 + i + 1) + ') = {' + str(i+1) + ', ' + str(200 + i + 2) + ', -' + str(100 + i + 1) + ', -' + str(200 + i + 1) +'};\n//+\n') #### poi mettere 201 al posto di 301

with open('wire.geo', 'a') as f:
    for i in range(N_POINTS):
        f.write('Plane Surface('+ str(2 + i + 1) + ') = {' + str(2 + i + 1) +'};\n//+\n')
"""