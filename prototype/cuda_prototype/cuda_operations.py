# -*- coding: utf-8 -*-

import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pycuda
import pycuda.driver as drv
import pycuda.autoinit
import numpy
from pycuda.compiler import SourceModule
from cuda_constant import *


# import pdb

def make_boundary_conditions_for_potentials_2(phi, n, boundary_type='50_50_BOUNDARY'):
    prev_phi = np.zeros([n[0], n[1], n[2]], dtype=np.float32)
    next_phi = np.zeros([n[0], n[1], n[2]], dtype=np.float32)
    ro = 1E-17
    if boundary_type == '50_50_BOUNDARY':
        for i in range(n[0]):
            for j in range(n[1]):
                for k in range(n[2]):
                        if i==n[0]-1 or j==n[1]-1 or k==n[2]-1 or i==0 or j==0 or k==0:
                            prev_phi[i][j][k] = phi
                            # next_phi[i][j][k] = 0
                            next_phi[i][j][k] = phi
                        else:
                            prev_phi[i][j][k] = 1#4*np.pi*ro*deltaT_
                            next_phi[i][j][k] = 0

    return (prev_phi, next_phi)


def potential_establish_method(prev_phi, next_phi, ro, epsilon=0.01):
    n = prev_phi.shape
    prev_sum, next_sum = epsilon, 3 * epsilon
    while np.abs(prev_sum - next_sum) > epsilon:
        prev_sum, next_sum = 0, 0
        for i in range(1, n[0] - 1):
            for j in range(1, n[1] - 1):
                for k in range(1, n[2] - 1):
                    next_phi[i][j][k] = deltaT_ * (
                        dif(prev_phi[i + 1][j][k], prev_phi[i][j][k], prev_phi[i - 1][j][k], xStepGrid) +
                        dif(prev_phi[i][j + 1][k], prev_phi[i][j][k], prev_phi[i][j - 1][k], yStepGrid) +
                        dif(prev_phi[i][j][k + 1], prev_phi[i][j][k], prev_phi[i][j][k - 1], zStepGrid) +
                        4 * np.pi * ro[i][j][k]) + prev_phi[i][j][k]
                    prev_sum += prev_phi[i][j][k]
                    next_sum += next_phi[i][j][k]
        for i in range(n[0]):
            for j in range(n[1]):
                for k in range(n[2]):
                    prev_phi[i][j][k] = next_phi[i][j][k]
        print('{}>{}'.format(np.abs(prev_sum - next_sum), epsilon))
    return (prev_phi, next_phi)


def potential_establish_method_2(prev_phi, next_phi, ro, epsilon=0.01):
    f = open('./kernel_code.c', 'r')
    mod = SourceModule("".join(f.readlines()))
    n = prev_phi.shape
    step = np.array([xStepGrid, yStepGrid, zStepGrid], dtype=np.float32)
    sizes = np.array(n, dtype=np.int32)
    potential_establish = mod.get_function("potential_establish")
    print(n[0],n[1],n[2])
    gd = int((n[0]-2)*(n[1]-2)*(n[2]-2)/100)
    if (gd != (n[0]-2)*(n[1]-2)*(n[2]-2)/100):
        gd += 1

    bd = int((n[0]-2)*(n[1]-2)*(n[2]-2)/gd)
    if (bd != (n[0]-2)*(n[1]-2)*(n[2]-2)/gd):
        bd += 1
    # gd, bd = 25, 390
    print('gd={}, bd = {} ===================='.format(gd, bd))

    prev_sum, next_sum = epsilon, 3 * epsilon
    while np.abs(prev_sum - next_sum) > epsilon:
        prev_sum, next_sum = 0, 0
        potential_establish(
            drv.InOut(prev_phi), drv.InOut(next_phi), drv.In(step), drv.In(sizes), 
            block=(bd,1, 1), grid=(gd,1))
        for i in range(1, n[0] - 1):
            for j in range(1, n[1] - 1):
                for k in range(1, n[2] - 1):
                    prev_sum += prev_phi[i][j][k]
                    next_sum += next_phi[i][j][k]
        # for kk in range(prev_phi.shape[0]):
        #     print(next_phi[kk])
            # plt.contourf(next_phi[kk])
            # plt.show()
        print('{}>{}'.format(np.abs(prev_sum - next_sum), epsilon))
        # for kk in range(next_phi.shape[0]):
        #     print(next_phi[kk])

    return (prev_phi, next_phi)

def potential_establish_method_3(prev_phi, next_phi, ro, epsilon=0.01):
    f = open('./kernel_code_one_block.c', 'r')
    mod = SourceModule("".join(f.readlines()))
    n = prev_phi.shape
    step = np.array([xStepGrid, yStepGrid, zStepGrid], dtype=np.float32)
    sizes = np.array(n, dtype=np.int32)
    potential_establish = mod.get_function("potential_establish")
    print(n[0],n[1],n[2])
    gd = int((n[0]-2)*(n[1]-2)*(n[2]-2)/512)
    if (gd != (n[0]-2)*(n[1]-2)*(n[2]-2)/512):
        gd += 1

    bd = int((n[0]-2)*(n[1]-2)*(n[2]-2)/gd)
    if (bd != (n[0]-2)*(n[1]-2)*(n[2]-2)/gd):
        bd += 1
    gd, bd = 1, 512
    print('gd={}, bd = {} ===================='.format(gd, bd))

    prev_sum, next_sum = epsilon, 3 * epsilon
    while np.abs(prev_sum - next_sum) > epsilon:
        prev_sum, next_sum = 0, 0
        potential_establish(
            drv.InOut(prev_phi), drv.InOut(next_phi), drv.In(step), drv.In(sizes), 
            block=(bd,1, 1), grid=(gd,1))
        for i in range(1, n[0] - 1):
            for j in range(1, n[1] - 1):
                for k in range(1, n[2] - 1):
                    prev_sum += prev_phi[i][j][k]
                    next_sum += next_phi[i][j][k]
        # for kk in range(prev_phi.shape[0]):
        #     print(next_phi[kk])
            # plt.contourf(next_phi[kk])
            # plt.show()
        print('{}>{}'.format(np.abs(prev_sum - next_sum), epsilon))
        # for kk in range(next_phi.shape[0]):
        #     print(next_phi[kk])

    return (prev_phi, next_phi)


def main():
    # Границы сетки
    xRange = np.linspace(0, xDimensionGrid, xNumberStepGrid + 1)
    yRange = np.linspace(0, yDimensionGrid, yNumberStepGrid + 1)
    zRange = np.linspace(0, zDimensionGrid, zNumberStepGrid + 1)
    # Граничные условия для потенциала
    n = (xRange.shape[0]+2, yRange.shape[0]+2, zRange.shape[0]+2)
    prev_phi, next_phi = make_boundary_conditions_for_potentials_2(50, n)
    for i in range(prev_phi.shape[0]):
        print(prev_phi[i])

    n = prev_phi.shape
    ro = np.zeros([n[0], n[1], n[2]], dtype=np.float32)
    # for i in range(1, n[0] - 1):
    #     for j in range(1, n[1] - 1):
    #         for k in range(1, n[2] - 1):
    #             ro[i][j][k] = chargeGridElectron[i-1][j-1][k-1] + chargeGridCarbon[i-1][j-1][k-1] + chargeGridHelium[i-1][j-1][k-1]
    print('THIS IS START OF potential_establish_method_2')
    prev_phi, next_phi = potential_establish_method_3(prev_phi, next_phi, ro, epsilon=0.01)
    print('THIS IS END OF potential_establish_method_2')

    for i in range(next_phi.shape[0]):
        # print(next_phi[i])
        cs = plt.contourf(next_phi[i])#, levels=[0, 10, 20, 30, 40, 45, 46, 47, 48, 49, 50, 51, 52])
        print(cs.levels)
        plt.colorbar(cs)
        plt.show()

if __name__ == '__main__':
    main()
