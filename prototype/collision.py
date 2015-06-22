# -*- coding: utf-8 -*-

import numpy as np
from time import time as timer

from constant import *

__author__ = 'fedul0x'

X_DIMENSION_COLLISION = 100
Y_DIMENSION_COLLISION = 100
Z_DIMENSION_COLLISION = 100

X_STEP_COLLISION = X_DIMENSION_GRID / X_DIMENSION_COLLISION
Y_STEP_COLLISION = Y_DIMENSION_GRID / Y_DIMENSION_COLLISION
Z_STEP_COLLISION = Z_DIMENSION_GRID / Z_DIMENSION_COLLISION

def make_collision_grid(x_dim, y_dim, z_dim):
    result = []
    for i in range(x_dim):
        result += [[]]
        for j in range(y_dim):
            result[i] += [[]]
            for k in range(z_dim):
                # result[i][j] += [[]]
                result[i][j] += [set()]
                
    return result

def distribute_by_grid(position, grid, time):
   # components = ((0, 0, 0), (-1, 0, 0), (0, -1, 0), (0, 0, -1),\
    #     (-1, -1, 0), (0, -1, -1), (-1, 0, -1), (-1, -1, -1),\
    #     (+1, 0, 0), (0, +1, 0), (0, 0, +1),\
    #     (+1, +1, 0), (0, +1, +1), (+1, 0, +1), (+1, +1, +1),)
    components = ((0, 0, 0), (-1, 0, 0), (0, -1, 0), (0, 0, -1),\
        (+1, 0, 0), (0, +1, 0), (0, 0, +1), )
        # (+1, +1, 0),     (-1, +1, 0),     (+1, -1, 0),     (-1, -1, 0), 
        # (+1, 0, +1),     (-1, 0, +1),     (+1, 0, -1),     (-1, 0, -1), 
        # (0, +1, +1),     (0, -1, +1),     (0, +1, -1),     (0, -1, -1),)
    n = position.shape[1]
    count = 0
    for num in range(n):
        x, y, z, radius, _, _, _, _ = position[time][num]
        for c in components:
            cx, cy, cz = x+c[0]*radius, y+c[1]*radius, z+c[2]*radius
            i, j, k = int(cx/X_STEP_COLLISION), int(cy/Y_STEP_COLLISION), int(cz/Z_STEP_COLLISION)
            if not (i<0 or j<0 or k<0 or i>X_DIMENSION_COLLISION-1 or j>Y_DIMENSION_COLLISION-1 or k>Z_DIMENSION_COLLISION-1):
                grid[i][j][k] = grid[i][j][k] | set([num])
                count += 1
            # n1, m1, o1 = int(cx/X_STEP_COLLISION), int(cy/Y_STEP_COLLISION), int(cz/Z_STEP_COLLISION)
            # n2, m2, o2 = int(x/X_STEP_COLLISION), int(y/Y_STEP_COLLISION), int(z/Z_STEP_COLLISION)
            # for i, j, k in zip(range(min(n1, n2), max(n1, n2)+1), range(min(m1, m2), max(m1, m2)+1), range(min(o1, o2), max(o1, o2)+1)):
            #     if not (i<0 or j<0 or k<0 or i>X_DIMENSION_COLLISION-1 or j>Y_DIMENSION_COLLISION-1 or k>Z_DIMENSION_COLLISION-1):
            #         grid[i][j][k] = grid[i][j][k] | set([num])
    return grid


def main():
    # Границы сетки
    x_range = np.linspace(0, X_DIMENSION_GRID, X_STEP_NUMBER_GRID + 1)
    y_range = np.linspace(0, Y_DIMENSION_GRID, Y_STEP_NUMBER_GRID + 1)
    z_range = np.linspace(0, Z_DIMENSION_GRID, Z_STEP_NUMBER_GRID + 1)
    n = (x_range.shape[0]-1) * (y_range.shape[0]-1) * (z_range.shape[0]-1)
    nc = int(n/(x_range.shape[0] - 1))
    # time, i, x, y, z, radius, charge, speedx, speedy, speedz
    electron = np.empty([DATA_IN_MEMORY_TIME, n, 6+2])
    carbon = \
        np.empty([DATA_IN_MEMORY_TIME, nc, 6+2])
    helium = np.empty([DATA_IN_MEMORY_TIME, n, 6+2])
    pe = electron
    ph = helium
    pc = carbon
    time = 0
    X_DIMENSION_GRID_2, Y_DIMENSION_GRID_2, Z_DIMENSION_GRID_2 = X_DIMENSION_GRID/10, Y_DIMENSION_GRID/10, Z_DIMENSION_GRID/10
    for i in range(n):
        electron[time][i][0] = np.random.rand() * X_DIMENSION_GRID_2
        electron[time][i][1] = np.random.rand() * Y_DIMENSION_GRID_2
        electron[time][i][2] = np.random.rand() * Z_DIMENSION_GRID_2
        electron[time][i][3] = 10E-6
        helium[time][i][0] = np.random.rand() * X_DIMENSION_GRID_2
        helium[time][i][1] = np.random.rand() * Y_DIMENSION_GRID_2
        helium[time][i][2] = np.random.rand() * Z_DIMENSION_GRID_2
        helium[time][i][3] = 10E-6
        if i <  nc:
            carbon[time][i][0] = 0
            carbon[time][i][1] = np.random.rand() * Y_DIMENSION_GRID_2
            carbon[time][i][2] = np.random.rand() * Z_DIMENSION_GRID_2
            carbon[time][i][3] = 10E-6

    count_1, count_2 = 0, 0

    start_time = timer()

    nc = carbon.shape[1]
    n = electron.shape[1]
    for i in range(nc):
        for j in range(n):
            xc, yc, zc, radiusc, _, _, _, _ = carbon[time][i] 
            xe, ye, ze, radiuse, _, _, _, _ = electron[time][j]
            # x, y, z, radius, _, _, _, _ = position[time][num]
            if np.sqrt((xc-xe)**2 + (yc-ye)**2 + (zc-ze)**2) < radiusc + radiuse:
                # print(i, j)
                count_1 += 1

    elapsed_time = timer() - start_time
    # print('===============')
    start_time_2 = timer()
    electron_collision = make_collision_grid(X_DIMENSION_COLLISION, Y_DIMENSION_COLLISION, Z_DIMENSION_COLLISION)
    helium_collision = make_collision_grid(X_DIMENSION_COLLISION, Y_DIMENSION_COLLISION, Z_DIMENSION_COLLISION)
    carbon_collision = make_collision_grid(X_DIMENSION_COLLISION, Y_DIMENSION_COLLISION, Z_DIMENSION_COLLISION)
    electron_collision = distribute_by_grid(electron, electron_collision, 0)
    helium_collision = distribute_by_grid(helium, helium_collision, 0)
    carbon_collision = distribute_by_grid(carbon, carbon_collision, 0)

    result = set()
    for i in range(X_DIMENSION_COLLISION):
        for j in range(Y_DIMENSION_COLLISION):
            for k in range(Z_DIMENSION_COLLISION):
                for a in carbon_collision[i][j][k]:
                    for b in electron_collision[i][j][k]:
                        xc, yc, zc, radiusc, _, _, _, _ = carbon[time][a] 
                        xe, ye, ze, radiuse, _, _, _, _ = electron[time][b]
                        # x, y, z, radius, _, _, _, _ = position[time][num]
                        if np.sqrt((xc-xe)**2 + (yc-ye)**2 + (zc-ze)**2) < radiusc + radiuse:
                            # print(a, b)
                            result |= set([(a, b,)])
                            count_2 += 1
    elapsed_time_2 = timer() - start_time_2

    # print(result)
    count_2 = len(result)
    print(elapsed_time)
    print(elapsed_time_2)
    print('{}={}'.format(count_1, count_2))
    # count = 0
    # for i in range(X_DIMENSION_COLLISION):
    #     for j in range(Y_DIMENSION_COLLISION):
    #         for k in range(Z_DIMENSION_COLLISION):
    #             count += len(electron_collision[i][j][k])
    # print(count)

    return (elapsed_time, elapsed_time_2)

if __name__ == '__main__':
    l1, l2 = [], []
    for i in range(5):
        a, b = main()
        l1 += [a]
        l2 += [b]
    print(l1)
    print(l2)

