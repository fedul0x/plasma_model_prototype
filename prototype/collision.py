# -*- coding: utf-8 -*-

import numpy as np
from time import time as timer
from math import *

from constant import *

__author__ = 'fedul0x'

# X_DIMENSION_COLLISION = 100
# Y_DIMENSION_COLLISION = 100
# Z_DIMENSION_COLLISION = 100

# X_STEP_COLLISION = X_DIMENSION_GRID / X_DIMENSION_COLLISION
# Y_STEP_COLLISION = Y_DIMENSION_GRID / Y_DIMENSION_COLLISION
# Z_STEP_COLLISION = Z_DIMENSION_GRID / Z_DIMENSION_COLLISION

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

def distribute_by_grid(position, grid, time, dims, steps):
    X_DIMENSION_COLLISION = dims[0]
    Y_DIMENSION_COLLISION = dims[1]
    Z_DIMENSION_COLLISION = dims[2]
    X_STEP_COLLISION = steps[0]
    Y_STEP_COLLISION = steps[1]
    Z_STEP_COLLISION = steps[2]


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

def check_collision(bp1, ep1, bp2, ep2, r1, r2, time_step):
    t0 = time_step
    x00, y00, z00 = bp1#0, 0, -10
    x01, y01, z01 = ep1#10, 10, 10
    x10, y10, z10 = bp2#0, 0, 10
    x11, y11, z11 = ep2#11, 10, -10
    # # r1, r2 = 1, 1
    # print('x00, y00, z00 = {}, {}, {}'.format(x00, y00, z00))
    # print('x01, y01, z01 = {}, {}, {}'.format(x01, y01, z01))
    # print('x10, y10, z10 = {}, {}, {}'.format(x10, y10, z10))
    # print('x11, y11, z11 = {}, {}, {}'.format(x11, y11, z11))
    # print('r1, r2 = {}, {}'.format(r1, r2))

    # # plt.plot([x00, x01], [y00, y01])
    # # plt.plot([x10, x11], [y10, y11])

    k0x = (x01 - x00) / t0
    k0y = (y01 - y00) / t0
    k0z = (z01 - z00) / t0
    k1z = (z11 - z10) / t0
    k1y = (y11 - y10) / t0
    k1x = (x11 - x10) / t0
    b0x, b0y, b0z = x00, y00, z00# - k0x*t0
    b1x, b1y, b1z = x10, y10, z10# - k0x*t0
    d = t0**2*(r1**2*x00**2 - 2*r1**2*x00*x01 - 2*r1**2*x00*x10 + 2*r1**2*x00*x11 + r1**2*x01**2 + 2*r1**2*x01*x10 - 2*r1**2*x01*x11 + r1**2*x10**2 - 2*r1**2*x10*x11 + r1**2*x11**2 + r1**2*y00**2 - 2*r1**2*y00*y01 - 2*r1**2*y00*y10 + 2*r1**2*y00*y11 + r1**2*y01**2 + 2*r1**2*y01*y10 - 2*r1**2*y01*y11 + r1**2*y10**2 - 2*r1**2*y10*y11 + r1**2*y11**2 + r1**2*z00**2 - 2*r1**2*z00*z01 - 2*r1**2*z00*z10 + 2*r1**2*z00*z11 + r1**2*z01**2 + 2*r1**2*z01*z10 - 2*r1**2*z01*z11 + r1**2*z10**2 - 2*r1**2*z10*z11 + r1**2*z11**2 + 2*r1*r2*x00**2 - 4*r1*r2*x00*x01 - 4*r1*r2*x00*x10 + 4*r1*r2*x00*x11 + 2*r1*r2*x01**2 + 4*r1*r2*x01*x10 - 4*r1*r2*x01*x11 + 2*r1*r2*x10**2 - 4*r1*r2*x10*x11 + 2*r1*r2*x11**2 + 2*r1*r2*y00**2 - 4*r1*r2*y00*y01 - 4*r1*r2*y00*y10 + 4*r1*r2*y00*y11 + 2*r1*r2*y01**2 + 4*r1*r2*y01*y10 - 4*r1*r2*y01*y11 + 2*r1*r2*y10**2 - 4*r1*r2*y10*y11 + 2*r1*r2*y11**2 + 2*r1*r2*z00**2 - 4*r1*r2*z00*z01 - 4*r1*r2*z00*z10 + 4*r1*r2*z00*z11 + 2*r1*r2*z01**2 + 4*r1*r2*z01*z10 - 4*r1*r2*z01*z11 + 2*r1*r2*z10**2 - 4*r1*r2*z10*z11 + 2*r1*r2*z11**2 + r2**2*x00**2 - 2*r2**2*x00*x01 - 2*r2**2*x00*x10 + 2*r2**2*x00*x11 + r2**2*x01**2 + 2*r2**2*x01*x10 - 2*r2**2*x01*x11 + r2**2*x10**2 - 2*r2**2*x10*x11 + r2**2*x11**2 + r2**2*y00**2 - 2*r2**2*y00*y01 - 2*r2**2*y00*y10 + 2*r2**2*y00*y11 + r2**2*y01**2 + 2*r2**2*y01*y10 - 2*r2**2*y01*y11 + r2**2*y10**2 - 2*r2**2*y10*y11 + r2**2*y11**2 + r2**2*z00**2 - 2*r2**2*z00*z01 - 2*r2**2*z00*z10 + 2*r2**2*z00*z11 + r2**2*z01**2 + 2*r2**2*z01*z10 - 2*r2**2*z01*z11 + r2**2*z10**2 - 2*r2**2*z10*z11 + r2**2*z11**2 - x00**2*y01**2 + 2*x00**2*y01*y11 - x00**2*y11**2 - x00**2*z01**2 + 2*x00**2*z01*z11 - x00**2*z11**2 + 2*x00*x01*y00*y01 - 2*x00*x01*y00*y11 - 2*x00*x01*y01*y10 + 2*x00*x01*y10*y11 + 2*x00*x01*z00*z01 - 2*x00*x01*z00*z11 - 2*x00*x01*z01*z10 + 2*x00*x01*z10*z11 + 2*x00*x10*y01**2 - 4*x00*x10*y01*y11 + 2*x00*x10*y11**2 + 2*x00*x10*z01**2 - 4*x00*x10*z01*z11 + 2*x00*x10*z11**2 - 2*x00*x11*y00*y01 + 2*x00*x11*y00*y11 + 2*x00*x11*y01*y10 - 2*x00*x11*y10*y11 - 2*x00*x11*z00*z01 + 2*x00*x11*z00*z11 + 2*x00*x11*z01*z10 - 2*x00*x11*z10*z11 - x01**2*y00**2 + 2*x01**2*y00*y10 - x01**2*y10**2 - x01**2*z00**2 + 2*x01**2*z00*z10 - x01**2*z10**2 - 2*x01*x10*y00*y01 + 2*x01*x10*y00*y11 + 2*x01*x10*y01*y10 - 2*x01*x10*y10*y11 - 2*x01*x10*z00*z01 + 2*x01*x10*z00*z11 + 2*x01*x10*z01*z10 - 2*x01*x10*z10*z11 + 2*x01*x11*y00**2 - 4*x01*x11*y00*y10 + 2*x01*x11*y10**2 + 2*x01*x11*z00**2 - 4*x01*x11*z00*z10 + 2*x01*x11*z10**2 - x10**2*y01**2 + 2*x10**2*y01*y11 - x10**2*y11**2 - x10**2*z01**2 + 2*x10**2*z01*z11 - x10**2*z11**2 + 2*x10*x11*y00*y01 - 2*x10*x11*y00*y11 - 2*x10*x11*y01*y10 + 2*x10*x11*y10*y11 + 2*x10*x11*z00*z01 - 2*x10*x11*z00*z11 - 2*x10*x11*z01*z10 + 2*x10*x11*z10*z11 - x11**2*y00**2 + 2*x11**2*y00*y10 - x11**2*y10**2 - x11**2*z00**2 + 2*x11**2*z00*z10 - x11**2*z10**2 - y00**2*z01**2 + 2*y00**2*z01*z11 - y00**2*z11**2 + 2*y00*y01*z00*z01 - 2*y00*y01*z00*z11 - 2*y00*y01*z01*z10 + 2*y00*y01*z10*z11 + 2*y00*y10*z01**2 - 4*y00*y10*z01*z11 + 2*y00*y10*z11**2 - 2*y00*y11*z00*z01 + 2*y00*y11*z00*z11 + 2*y00*y11*z01*z10 - 2*y00*y11*z10*z11 - y01**2*z00**2 + 2*y01**2*z00*z10 - y01**2*z10**2 - 2*y01*y10*z00*z01 + 2*y01*y10*z00*z11 + 2*y01*y10*z01*z10 - 2*y01*y10*z10*z11 + 2*y01*y11*z00**2 - 4*y01*y11*z00*z10 + 2*y01*y11*z10**2 - y10**2*z01**2 + 2*y10**2*z01*z11 - y10**2*z11**2 + 2*y10*y11*z00*z01 - 2*y10*y11*z00*z11 - 2*y10*y11*z01*z10 + 2*y10*y11*z10*z11 - y11**2*z00**2 + 2*y11**2*z00*z10 - y11**2*z10**2)
    if d >= 0:
        t = (t0*x00**2 - t0*x00*x01 - 2*t0*x00*x10 + t0*x00*x11 + t0*x01*x10 + t0*x10**2 - t0*x10*x11 + t0*y00**2 - t0*y00*y01 - 2*t0*y00*y10 + t0*y00*y11 + t0*y01*y10 + t0*y10**2 - t0*y10*y11 + t0*z00**2 - t0*z00*z01 - 2*t0*z00*z10 + t0*z00*z11 + t0*z01*z10 + t0*z10**2 - t0*z10*z11 - sqrt(d))/(x00**2 - 2*x00*x01 - 2*x00*x10 + 2*x00*x11 + x01**2 + 2*x01*x10 - 2*x01*x11 + x10**2 - 2*x10*x11 + x11**2 + y00**2 - 2*y00*y01 - 2*y00*y10 + 2*y00*y11 + y01**2 + 2*y01*y10 - 2*y01*y11 + y10**2 - 2*y10*y11 + y11**2 + z00**2 - 2*z00*z01 - 2*z00*z10 + 2*z00*z11 + z01**2 + 2*z01*z10 - 2*z01*z11 + z10**2 - 2*z10*z11 + z11**2)
        print(d, t0)
        if t >= 0 and t <= t0:
            return (t0, None)
    
    return (None, None)

    # v = (t0*(k1x*t0*x00 - k1x*t0*x10 + k1y*t0*y00 - k1y*t0*y10 + k1z*t0*z00 - k1z*t0*z10 + x00**2 - x00*x01 - x00*x10 + x01*x10 + y00**2 - y00*y01 - y00*y10 + y01*y10 + z00**2 - z00*z01 - z00*z10 + z01*z10) + sqrt(t0**2*(k1x**2*r1**2*t0**2 + 2*k1x**2*r1*r2*t0**2 + k1x**2*r2**2*t0**2 - k1x**2*t0**2*y00**2 + 2*k1x**2*t0**2*y00*y10 - k1x**2*t0**2*y10**2 - k1x**2*t0**2*z00**2 + 2*k1x**2*t0**2*z00*z10 - k1x**2*t0**2*z10**2 + 2*k1x*k1y*t0**2*x00*y00 - 2*k1x*k1y*t0**2*x00*y10 - 2*k1x*k1y*t0**2*x10*y00 + 2*k1x*k1y*t0**2*x10*y10 + 2*k1x*k1z*t0**2*x00*z00 - 2*k1x*k1z*t0**2*x00*z10 - 2*k1x*k1z*t0**2*x10*z00 + 2*k1x*k1z*t0**2*x10*z10 + 2*k1x*r1**2*t0*x00 - 2*k1x*r1**2*t0*x01 + 4*k1x*r1*r2*t0*x00 - 4*k1x*r1*r2*t0*x01 + 2*k1x*r2**2*t0*x00 - 2*k1x*r2**2*t0*x01 - 2*k1x*t0*x00*y00*y01 + 2*k1x*t0*x00*y00*y10 + 2*k1x*t0*x00*y01*y10 - 2*k1x*t0*x00*y10**2 - 2*k1x*t0*x00*z00*z01 + 2*k1x*t0*x00*z00*z10 + 2*k1x*t0*x00*z01*z10 - 2*k1x*t0*x00*z10**2 + 2*k1x*t0*x01*y00**2 - 4*k1x*t0*x01*y00*y10 + 2*k1x*t0*x01*y10**2 + 2*k1x*t0*x01*z00**2 - 4*k1x*t0*x01*z00*z10 + 2*k1x*t0*x01*z10**2 - 2*k1x*t0*x10*y00**2 + 2*k1x*t0*x10*y00*y01 + 2*k1x*t0*x10*y00*y10 - 2*k1x*t0*x10*y01*y10 - 2*k1x*t0*x10*z00**2 + 2*k1x*t0*x10*z00*z01 + 2*k1x*t0*x10*z00*z10 - 2*k1x*t0*x10*z01*z10 + k1y**2*r1**2*t0**2 + 2*k1y**2*r1*r2*t0**2 + k1y**2*r2**2*t0**2 - k1y**2*t0**2*x00**2 + 2*k1y**2*t0**2*x00*x10 - k1y**2*t0**2*x10**2 - k1y**2*t0**2*z00**2 + 2*k1y**2*t0**2*z00*z10 - k1y**2*t0**2*z10**2 + 2*k1y*k1z*t0**2*y00*z00 - 2*k1y*k1z*t0**2*y00*z10 - 2*k1y*k1z*t0**2*y10*z00 + 2*k1y*k1z*t0**2*y10*z10 + 2*k1y*r1**2*t0*y00 - 2*k1y*r1**2*t0*y01 + 4*k1y*r1*r2*t0*y00 - 4*k1y*r1*r2*t0*y01 + 2*k1y*r2**2*t0*y00 - 2*k1y*r2**2*t0*y01 + 2*k1y*t0*x00**2*y01 - 2*k1y*t0*x00**2*y10 - 2*k1y*t0*x00*x01*y00 + 2*k1y*t0*x00*x01*y10 + 2*k1y*t0*x00*x10*y00 - 4*k1y*t0*x00*x10*y01 + 2*k1y*t0*x00*x10*y10 + 2*k1y*t0*x01*x10*y00 - 2*k1y*t0*x01*x10*y10 - 2*k1y*t0*x10**2*y00 + 2*k1y*t0*x10**2*y01 - 2*k1y*t0*y00*z00*z01 + 2*k1y*t0*y00*z00*z10 + 2*k1y*t0*y00*z01*z10 - 2*k1y*t0*y00*z10**2 + 2*k1y*t0*y01*z00**2 - 4*k1y*t0*y01*z00*z10 + 2*k1y*t0*y01*z10**2 - 2*k1y*t0*y10*z00**2 + 2*k1y*t0*y10*z00*z01 + 2*k1y*t0*y10*z00*z10 - 2*k1y*t0*y10*z01*z10 + k1z**2*r1**2*t0**2 + 2*k1z**2*r1*r2*t0**2 + k1z**2*r2**2*t0**2 - k1z**2*t0**2*x00**2 + 2*k1z**2*t0**2*x00*x10 - k1z**2*t0**2*x10**2 - k1z**2*t0**2*y00**2 + 2*k1z**2*t0**2*y00*y10 - k1z**2*t0**2*y10**2 + 2*k1z*r1**2*t0*z00 - 2*k1z*r1**2*t0*z01 + 4*k1z*r1*r2*t0*z00 - 4*k1z*r1*r2*t0*z01 + 2*k1z*r2**2*t0*z00 - 2*k1z*r2**2*t0*z01 + 2*k1z*t0*x00**2*z01 - 2*k1z*t0*x00**2*z10 - 2*k1z*t0*x00*x01*z00 + 2*k1z*t0*x00*x01*z10 + 2*k1z*t0*x00*x10*z00 - 4*k1z*t0*x00*x10*z01 + 2*k1z*t0*x00*x10*z10 + 2*k1z*t0*x01*x10*z00 - 2*k1z*t0*x01*x10*z10 - 2*k1z*t0*x10**2*z00 + 2*k1z*t0*x10**2*z01 + 2*k1z*t0*y00**2*z01 - 2*k1z*t0*y00**2*z10 - 2*k1z*t0*y00*y01*z00 + 2*k1z*t0*y00*y01*z10 + 2*k1z*t0*y00*y10*z00 - 4*k1z*t0*y00*y10*z01 + 2*k1z*t0*y00*y10*z10 + 2*k1z*t0*y01*y10*z00 - 2*k1z*t0*y01*y10*z10 - 2*k1z*t0*y10**2*z00 + 2*k1z*t0*y10**2*z01 + r1**2*x00**2 - 2*r1**2*x00*x01 + r1**2*x01**2 + r1**2*y00**2 - 2*r1**2*y00*y01 + r1**2*y01**2 + r1**2*z00**2 - 2*r1**2*z00*z01 + r1**2*z01**2 + 2*r1*r2*x00**2 - 4*r1*r2*x00*x01 + 2*r1*r2*x01**2 + 2*r1*r2*y00**2 - 4*r1*r2*y00*y01 + 2*r1*r2*y01**2 + 2*r1*r2*z00**2 - 4*r1*r2*z00*z01 + 2*r1*r2*z01**2 + r2**2*x00**2 - 2*r2**2*x00*x01 + r2**2*x01**2 + r2**2*y00**2 - 2*r2**2*y00*y01 + r2**2*y01**2 + r2**2*z00**2 - 2*r2**2*z00*z01 + r2**2*z01**2 - x00**2*y01**2 + 2*x00**2*y01*y10 - x00**2*y10**2 - x00**2*z01**2 + 2*x00**2*z01*z10 - x00**2*z10**2 + 2*x00*x01*y00*y01 - 2*x00*x01*y00*y10 - 2*x00*x01*y01*y10 + 2*x00*x01*y10**2 + 2*x00*x01*z00*z01 - 2*x00*x01*z00*z10 - 2*x00*x01*z01*z10 + 2*x00*x01*z10**2 - 2*x00*x10*y00*y01 + 2*x00*x10*y00*y10 + 2*x00*x10*y01**2 - 2*x00*x10*y01*y10 - 2*x00*x10*z00*z01 + 2*x00*x10*z00*z10 + 2*x00*x10*z01**2 - 2*x00*x10*z01*z10 - x01**2*y00**2 + 2*x01**2*y00*y10 - x01**2*y10**2 - x01**2*z00**2 + 2*x01**2*z00*z10 - x01**2*z10**2 + 2*x01*x10*y00**2 - 2*x01*x10*y00*y01 - 2*x01*x10*y00*y10 + 2*x01*x10*y01*y10 + 2*x01*x10*z00**2 - 2*x01*x10*z00*z01 - 2*x01*x10*z00*z10 + 2*x01*x10*z01*z10 - x10**2*y00**2 + 2*x10**2*y00*y01 - x10**2*y01**2 - x10**2*z00**2 + 2*x10**2*z00*z01 - x10**2*z01**2 - y00**2*z01**2 + 2*y00**2*z01*z10 - y00**2*z10**2 + 2*y00*y01*z00*z01 - 2*y00*y01*z00*z10 - 2*y00*y01*z01*z10 + 2*y00*y01*z10**2 - 2*y00*y10*z00*z01 + 2*y00*y10*z00*z10 + 2*y00*y10*z01**2 - 2*y00*y10*z01*z10 - y01**2*z00**2 + 2*y01**2*z00*z10 - y01**2*z10**2 + 2*y01*y10*z00**2 - 2*y01*y10*z00*z01 - 2*y01*y10*z00*z10 + 2*y01*y10*z01*z10 - y10**2*z00**2 + 2*y10**2*z00*z01 - y10**2*z01**2)))/(k1x**2*t0**2 + 2*k1x*t0*x00 - 2*k1x*t0*x01 + k1y**2*t0**2 + 2*k1y*t0*y00 - 2*k1y*t0*y01 + k1z**2*t0**2 + 2*k1z*t0*z00 - 2*k1z*t0*z01 + x00**2 - 2*x00*x01 + x01**2 + y00**2 - 2*y00*y01 + y01**2 + z00**2 - 2*z00*z01 + z01**2)
    # d = t0**2*(r1**2*x00**2 - 2*r1**2*x00*x01 - 2*r1**2*x00*x10 + 2*r1**2*x00*x11 + r1**2*x01**2 + 2*r1**2*x01*x10 - 2*r1**2*x01*x11 + r1**2*x10**2 - 2*r1**2*x10*x11 + r1**2*x11**2 + r1**2*y00**2 - 2*r1**2*y00*y01 - 2*r1**2*y00*y10 + 2*r1**2*y00*y11 + r1**2*y01**2 + 2*r1**2*y01*y10 - 2*r1**2*y01*y11 + r1**2*y10**2 - 2*r1**2*y10*y11 + r1**2*y11**2 + r1**2*z00**2 - 2*r1**2*z00*z01 - 2*r1**2*z00*z10 + 2*r1**2*z00*z11 + r1**2*z01**2 + 2*r1**2*z01*z10 - 2*r1**2*z01*z11 + r1**2*z10**2 - 2*r1**2*z10*z11 + r1**2*z11**2 + 2*r1*r2*x00**2 - 4*r1*r2*x00*x01 - 4*r1*r2*x00*x10 + 4*r1*r2*x00*x11 + 2*r1*r2*x01**2 + 4*r1*r2*x01*x10 - 4*r1*r2*x01*x11 + 2*r1*r2*x10**2 - 4*r1*r2*x10*x11 + 2*r1*r2*x11**2 + 2*r1*r2*y00**2 - 4*r1*r2*y00*y01 - 4*r1*r2*y00*y10 + 4*r1*r2*y00*y11 + 2*r1*r2*y01**2 + 4*r1*r2*y01*y10 - 4*r1*r2*y01*y11 + 2*r1*r2*y10**2 - 4*r1*r2*y10*y11 + 2*r1*r2*y11**2 + 2*r1*r2*z00**2 - 4*r1*r2*z00*z01 - 4*r1*r2*z00*z10 + 4*r1*r2*z00*z11 + 2*r1*r2*z01**2 + 4*r1*r2*z01*z10 - 4*r1*r2*z01*z11 + 2*r1*r2*z10**2 - 4*r1*r2*z10*z11 + 2*r1*r2*z11**2 + r2**2*x00**2 - 2*r2**2*x00*x01 - 2*r2**2*x00*x10 + 2*r2**2*x00*x11 + r2**2*x01**2 + 2*r2**2*x01*x10 - 2*r2**2*x01*x11 + r2**2*x10**2 - 2*r2**2*x10*x11 + r2**2*x11**2 + r2**2*y00**2 - 2*r2**2*y00*y01 - 2*r2**2*y00*y10 + 2*r2**2*y00*y11 + r2**2*y01**2 + 2*r2**2*y01*y10 - 2*r2**2*y01*y11 + r2**2*y10**2 - 2*r2**2*y10*y11 + r2**2*y11**2 + r2**2*z00**2 - 2*r2**2*z00*z01 - 2*r2**2*z00*z10 + 2*r2**2*z00*z11 + r2**2*z01**2 + 2*r2**2*z01*z10 - 2*r2**2*z01*z11 + r2**2*z10**2 - 2*r2**2*z10*z11 + r2**2*z11**2 - x00**2*y01**2 + 2*x00**2*y01*y11 - x00**2*y11**2 - x00**2*z01**2 + 2*x00**2*z01*z11 - x00**2*z11**2 + 2*x00*x01*y00*y01 - 2*x00*x01*y00*y11 - 2*x00*x01*y01*y10 + 2*x00*x01*y10*y11 + 2*x00*x01*z00*z01 - 2*x00*x01*z00*z11 - 2*x00*x01*z01*z10 + 2*x00*x01*z10*z11 + 2*x00*x10*y01**2 - 4*x00*x10*y01*y11 + 2*x00*x10*y11**2 + 2*x00*x10*z01**2 - 4*x00*x10*z01*z11 + 2*x00*x10*z11**2 - 2*x00*x11*y00*y01 + 2*x00*x11*y00*y11 + 2*x00*x11*y01*y10 - 2*x00*x11*y10*y11 - 2*x00*x11*z00*z01 + 2*x00*x11*z00*z11 + 2*x00*x11*z01*z10 - 2*x00*x11*z10*z11 - x01**2*y00**2 + 2*x01**2*y00*y10 - x01**2*y10**2 - x01**2*z00**2 + 2*x01**2*z00*z10 - x01**2*z10**2 - 2*x01*x10*y00*y01 + 2*x01*x10*y00*y11 + 2*x01*x10*y01*y10 - 2*x01*x10*y10*y11 - 2*x01*x10*z00*z01 + 2*x01*x10*z00*z11 + 2*x01*x10*z01*z10 - 2*x01*x10*z10*z11 + 2*x01*x11*y00**2 - 4*x01*x11*y00*y10 + 2*x01*x11*y10**2 + 2*x01*x11*z00**2 - 4*x01*x11*z00*z10 + 2*x01*x11*z10**2 - x10**2*y01**2 + 2*x10**2*y01*y11 - x10**2*y11**2 - x10**2*z01**2 + 2*x10**2*z01*z11 - x10**2*z11**2 + 2*x10*x11*y00*y01 - 2*x10*x11*y00*y11 - 2*x10*x11*y01*y10 + 2*x10*x11*y10*y11 + 2*x10*x11*z00*z01 - 2*x10*x11*z00*z11 - 2*x10*x11*z01*z10 + 2*x10*x11*z10*z11 - x11**2*y00**2 + 2*x11**2*y00*y10 - x11**2*y10**2 - x11**2*z00**2 + 2*x11**2*z00*z10 - x11**2*z10**2 - y00**2*z01**2 + 2*y00**2*z01*z11 - y00**2*z11**2 + 2*y00*y01*z00*z01 - 2*y00*y01*z00*z11 - 2*y00*y01*z01*z10 + 2*y00*y01*z10*z11 + 2*y00*y10*z01**2 - 4*y00*y10*z01*z11 + 2*y00*y10*z11**2 - 2*y00*y11*z00*z01 + 2*y00*y11*z00*z11 + 2*y00*y11*z01*z10 - 2*y00*y11*z10*z11 - y01**2*z00**2 + 2*y01**2*z00*z10 - y01**2*z10**2 - 2*y01*y10*z00*z01 + 2*y01*y10*z00*z11 + 2*y01*y10*z01*z10 - 2*y01*y10*z10*z11 + 2*y01*y11*z00**2 - 4*y01*y11*z00*z10 + 2*y01*y11*z10**2 - y10**2*z01**2 + 2*y10**2*z01*z11 - y10**2*z11**2 + 2*y10*y11*z00*z01 - 2*y10*y11*z00*z11 - 2*y10*y11*z01*z10 + 2*y10*y11*z10*z11 - y11**2*z00**2 + 2*y11**2*z00*z10 - y11**2*z10**2)
    # d1 = t0**2*(r1**2*x00**2 - 2*r1**2*x00*x01 - 2*r1**2*x00*x10 + 2*r1**2*x00*x11 + r1**2*x01**2 + 2*r1**2*x01*x10 - 2*r1**2*x01*x11 + r1**2*x10**2 - 2*r1**2*x10*x11 + r1**2*x11**2 + r1**2*y00**2 - 2*r1**2*y00*y01 - 2*r1**2*y00*y10 + 2*r1**2*y00*y11 + r1**2*y01**2 + 2*r1**2*y01*y10 - 2*r1**2*y01*y11 + r1**2*y10**2 - 2*r1**2*y10*y11 + r1**2*y11**2 + r1**2*z00**2 - 2*r1**2*z00*z01 - 2*r1**2*z00*z10 + 2*r1**2*z00*z11 + r1**2*z01**2 + 2*r1**2*z01*z10 - 2*r1**2*z01*z11 + r1**2*z10**2 - 2*r1**2*z10*z11 + r1**2*z11**2 + 2*r1*r2*x00**2 - 4*r1*r2*x00*x01 - 4*r1*r2*x00*x10 + 4*r1*r2*x00*x11 + 2*r1*r2*x01**2 + 4*r1*r2*x01*x10 - 4*r1*r2*x01*x11 + 2*r1*r2*x10**2 - 4*r1*r2*x10*x11 + 2*r1*r2*x11**2 + 2*r1*r2*y00**2 - 4*r1*r2*y00*y01 - 4*r1*r2*y00*y10 + 4*r1*r2*y00*y11 + 2*r1*r2*y01**2 + 4*r1*r2*y01*y10 - 4*r1*r2*y01*y11 + 2*r1*r2*y10**2 - 4*r1*r2*y10*y11 + 2*r1*r2*y11**2 + 2*r1*r2*z00**2 - 4*r1*r2*z00*z01 - 4*r1*r2*z00*z10 + 4*r1*r2*z00*z11 + 2*r1*r2*z01**2 + 4*r1*r2*z01*z10 - 4*r1*r2*z01*z11 + 2*r1*r2*z10**2 - 4*r1*r2*z10*z11 + 2*r1*r2*z11**2 + r2**2*x00**2 - 2*r2**2*x00*x01 - 2*r2**2*x00*x10 + 2*r2**2*x00*x11 + r2**2*x01**2 + 2*r2**2*x01*x10 - 2*r2**2*x01*x11 + r2**2*x10**2 - 2*r2**2*x10*x11 + r2**2*x11**2 + r2**2*y00**2 - 2*r2**2*y00*y01 - 2*r2**2*y00*y10 + 2*r2**2*y00*y11 + r2**2*y01**2 + 2*r2**2*y01*y10 - 2*r2**2*y01*y11 + r2**2*y10**2 - 2*r2**2*y10*y11 + r2**2*y11**2 + r2**2*z00**2 - 2*r2**2*z00*z01 - 2*r2**2*z00*z10 + 2*r2**2*z00*z11 + r2**2*z01**2 + 2*r2**2*z01*z10 - 2*r2**2*z01*z11 + r2**2*z10**2 - 2*r2**2*z10*z11 + r2**2*z11**2 - x00**2*y01**2 + 2*x00**2*y01*y11 - x00**2*y11**2 - x00**2*z01**2 + 2*x00**2*z01*z11 - x00**2*z11**2 + 2*x00*x01*y00*y01 - 2*x00*x01*y00*y11 - 2*x00*x01*y01*y10 + 2*x00*x01*y10*y11 + 2*x00*x01*z00*z01 - 2*x00*x01*z00*z11 - 2*x00*x01*z01*z10 + 2*x00*x01*z10*z11 + 2*x00*x10*y01**2 - 4*x00*x10*y01*y11 + 2*x00*x10*y11**2 + 2*x00*x10*z01**2 - 4*x00*x10*z01*z11 + 2*x00*x10*z11**2 - 2*x00*x11*y00*y01 + 2*x00*x11*y00*y11 + 2*x00*x11*y01*y10 - 2*x00*x11*y10*y11 - 2*x00*x11*z00*z01 + 2*x00*x11*z00*z11 + 2*x00*x11*z01*z10 - 2*x00*x11*z10*z11 - x01**2*y00**2 + 2*x01**2*y00*y10 - x01**2*y10**2 - x01**2*z00**2 + 2*x01**2*z00*z10 - x01**2*z10**2 - 2*x01*x10*y00*y01 + 2*x01*x10*y00*y11 + 2*x01*x10*y01*y10 - 2*x01*x10*y10*y11 - 2*x01*x10*z00*z01 + 2*x01*x10*z00*z11 + 2*x01*x10*z01*z10 - 2*x01*x10*z10*z11 + 2*x01*x11*y00**2 - 4*x01*x11*y00*y10 + 2*x01*x11*y10**2 + 2*x01*x11*z00**2 - 4*x01*x11*z00*z10 + 2*x01*x11*z10**2 - x10**2*y01**2 + 2*x10**2*y01*y11 - x10**2*y11**2 - x10**2*z01**2 + 2*x10**2*z01*z11 - x10**2*z11**2 + 2*x10*x11*y00*y01 - 2*x10*x11*y00*y11 - 2*x10*x11*y01*y10 + 2*x10*x11*y10*y11 + 2*x10*x11*z00*z01 - 2*x10*x11*z00*z11 - 2*x10*x11*z01*z10 + 2*x10*x11*z10*z11 - x11**2*y00**2 + 2*x11**2*y00*y10 - x11**2*y10**2 - x11**2*z00**2 + 2*x11**2*z00*z10 - x11**2*z10**2 - y00**2*z01**2 + 2*y00**2*z01*z11 - y00**2*z11**2 + 2*y00*y01*z00*z01 - 2*y00*y01*z00*z11 - 2*y00*y01*z01*z10 + 2*y00*y01*z10*z11 + 2*y00*y10*z01**2 - 4*y00*y10*z01*z11 + 2*y00*y10*z11**2 - 2*y00*y11*z00*z01 + 2*y00*y11*z00*z11 + 2*y00*y11*z01*z10 - 2*y00*y11*z10*z11 - y01**2*z00**2 + 2*y01**2*z00*z10 - y01**2*z10**2 - 2*y01*y10*z00*z01 + 2*y01*y10*z00*z11 + 2*y01*y10*z01*z10 - 2*y01*y10*z10*z11 + 2*y01*y11*z00**2 - 4*y01*y11*z00*z10 + 2*y01*y11*z10**2 - y10**2*z01**2 + 2*y10**2*z01*z11 - y10**2*z11**2 + 2*y10*y11*z00*z01 - 2*y10*y11*z00*z11 - 2*y10*y11*z01*z10 + 2*y10*y11*z10*z11 - y11**2*z00**2 + 2*y11**2*z00*z10 - y11**2*z10**2)
    # assert(d==d1)
    # print('=========================\n{}\n{}\n========================='.format(d, d1))
    # try:
    #     v_1 = (t0*x00**2 - t0*x00*x01 - 2*t0*x00*x10 + t0*x00*x11 + t0*x01*x10 + t0*x10**2 - t0*x10*x11 + t0*y00**2 - t0*y00*y01 - 2*t0*y00*y10 + t0*y00*y11 + t0*y01*y10 + t0*y10**2 - t0*y10*y11 + t0*z00**2 - t0*z00*z01 - 2*t0*z00*z10 + t0*z00*z11 + t0*z01*z10 + t0*z10**2 - t0*z10*z11 - sqrt(d))/(x00**2 - 2*x00*x01 - 2*x00*x10 + 2*x00*x11 + x01**2 + 2*x01*x10 - 2*x01*x11 + x10**2 - 2*x10*x11 + x11**2 + y00**2 - 2*y00*y01 - 2*y00*y10 + 2*y00*y11 + y01**2 + 2*y01*y10 - 2*y01*y11 + y10**2 - 2*y10*y11 + y11**2 + z00**2 - 2*z00*z01 - 2*z00*z10 + 2*z00*z11 + z01**2 + 2*z01*z10 - 2*z01*z11 + z10**2 - 2*z10*z11 + z11**2)
    # except ValueError:
    #     v_1 = None
    # try:
    #     v_2 = (t0*(x00**2 - x00*x01 - 2*x00*x10 + x00*x11 + x01*x10 + x10**2 - x10*x11 + y00**2 - y00*y01 - 2*y00*y10 + y00*y11 + y01*y10 + y10**2 - y10*y11 + z00**2 - z00*z01 - 2*z00*z10 + z00*z11 + z01*z10 + z10**2 - z10*z11) + sqrt(d1))/(x00**2 - 2*x00*x01 - 2*x00*x10 + 2*x00*x11 + x01**2 + 2*x01*x10 - 2*x01*x11 + x10**2 - 2*x10*x11 + x11**2 + y00**2 - 2*y00*y01 - 2*y00*y10 + 2*y00*y11 + y01**2 + 2*y01*y10 - 2*y01*y11 + y10**2 - 2*y10*y11 + y11**2 + z00**2 - 2*z00*z01 - 2*z00*z10 + 2*z00*z11 + z01**2 + 2*z01*z10 - 2*z01*z11 + z10**2 - 2*z10*z11 + z11**2)
    # except ValueError:
    #     v_2 = None
    # # print('CRASH 1:\n[{}, {}]'.format(v_1, v_2))
    # vs_2 = [None, None]
    # if not (v_1 == None and v_2 == None):
    #     vs = sorted([v_1, v_2])
    #     # print('CRASH 2:\n{}'.format(vs))
    #     for i, j in zip(vs, range(len(vs_2))):
    #         if i >=0 and i<=t0:
    #             vs_2[j] = i
    #     # print('CRASH 3:\n{}'.format(vs_2))

    # return vs_2


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

