import numpy as np

N_POINTS = 1500
R_EXT = 15

with open('spira.geo', 'w') as f:
    for i in range(N_POINTS):
        f.write('Point('+ str(i+1) + ') = {' + str(R_EXT*np.cos(2*np.pi/N_POINTS*i)) + ', ' + str(R_EXT*np.sin(2*np.pi/N_POINTS*i)) + ', 0, 1.0};\n//+\n')

with open('spira.geo', 'a') as f:
    for i in range(N_POINTS-1):
        f.write('Line('+ str(i+1) + ') = {' + str(i+1) + ', ' + str(i+2) +'};\n//+\n')
    f.write('Line(' + str(N_POINTS) + ') = {' + str(N_POINTS) + ', 1};\n//+\n')

with open('spira.geo', 'a') as f:
    f.write('Physical Curve("spira") = {' + str([i for i in range(1, N_POINTS+1)]) + '};\n//+\n')

"""
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100};
//+
Plane Surface(1) = {1};
//+
Physical Surface("domain") = {1};
"""