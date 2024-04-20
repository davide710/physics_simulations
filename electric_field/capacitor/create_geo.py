"""
######################## faccia sotto
base_x = 0
base_y = 0
base_z = 0

length = 10
width = 10
height = 0.1

vertices = [
    (base_x, base_y, base_z),
    (base_x + length, base_y, base_z),
    (base_x + length, base_y + width, base_z),
    (base_x, base_y + width, base_z),
    (base_x, base_y, base_z + height),
    (base_x + length, base_y, base_z + height),
    (base_x + length, base_y + width, base_z + height),
    (base_x, base_y + width, base_z + height),
]

with open("capacitor.geo", "w") as file:
    for i, vertex in enumerate(vertices):
        file.write(f'Point({i+1}) = ' + '{' + str(vertex[0]) + ', ' + str(vertex[1]) + ', ' + str(vertex[2]) + '};\n//+\n')

    for i in range(4):
        file.write(f'Line({i+1}) = ' + '{' + str(i+1) + ', ' + str((i+1)%4 + 1) + '};\n//+\n')
        file.write(f'Line({i + 5}) = ' + '{' + str(i+5) + ', ' + str((i+1)%4 + 5) + '};\n//+\n')

    for i in range(4):
        file.write(f'Line({i+9}) = ' + '{' + str(i+1) + ', ' + str(i+5) + '};\n//+\n')
"""

######################## faccia sopra
base_x = 0
base_y = 0
base_z = 0.3

length = 10
width = 10
height = 0.1

vertices = [
    (base_x, base_y, base_z),
    (base_x + length, base_y, base_z),
    (base_x + length, base_y + width, base_z),
    (base_x, base_y + width, base_z),
    (base_x, base_y, base_z + height),
    (base_x + length, base_y, base_z + height),
    (base_x + length, base_y + width, base_z + height),
    (base_x, base_y + width, base_z + height),
]

with open("capacitor.geo", "a") as file:
    for i, vertex in enumerate(vertices):
        file.write(f'Point({i+9}) = ' + '{' + str(vertex[0]) + ', ' + str(vertex[1]) + ', ' + str(vertex[2]) + '};\n//+\n')

    for i in range(4):
        file.write(f'Line({i+13}) = ' + '{' + str(i+1 +8) + ', ' + str((i+1)%4 + 1 +8) + '};\n//+\n')
        file.write(f'Line({i + 17}) = ' + '{' + str(i+5 +8) + ', ' + str((i+1)%4 + 5 +8) + '};\n//+\n')

    for i in range(4):
        file.write(f'Line({i+21}) = ' + '{' + str(i+1 +8) + ', ' + str(i+5 +8) + '};\n//+\n')
