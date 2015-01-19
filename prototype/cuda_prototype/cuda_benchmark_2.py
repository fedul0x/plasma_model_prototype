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
import time

def dif(x, y, z, step):
    return (x - 2 * y + z) / step ** 2

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
    iter_count = 0
    ro = 1E-17
    while np.abs(prev_sum - next_sum) > epsilon:
        prev_sum, next_sum = 0, 0
        for i in range(1, n[0] - 1):
            for j in range(1, n[1] - 1):
                for k in range(1, n[2] - 1):
                    next_phi[i][j][k] = deltaT_ * (
                        dif(prev_phi[i + 1][j][k], prev_phi[i][j][k], prev_phi[i - 1][j][k], xStepGrid) +
                        dif(prev_phi[i][j + 1][k], prev_phi[i][j][k], prev_phi[i][j - 1][k], yStepGrid) +
                        dif(prev_phi[i][j][k + 1], prev_phi[i][j][k], prev_phi[i][j][k - 1], zStepGrid) +
                        4 * np.pi * ro) + prev_phi[i][j][k]
                    prev_sum += prev_phi[i][j][k]
                    next_sum += next_phi[i][j][k]
        prev_phi, next_phi = next_phi, prev_phi
        # for i in range(n[0]):
        #     for j in range(n[1]):
        #         for k in range(n[2]):
        #             prev_phi[i][j][k] = next_phi[i][j][k]
        iter_count += 1
        # print('{}>{}'.format(np.abs(prev_sum - next_sum), epsilon))
    print('Iteration count = {}'.format(iter_count))

    return (prev_phi, next_phi)

def potential_establish_method_cuda(file_name, prev_phi, next_phi, ro, epsilon=0.01):
    f = open(file_name, 'r')
    mod = SourceModule("".join(f.readlines()))#, options=['-ptx'])
    n = prev_phi.shape
    step = np.array([xStepGrid, yStepGrid, zStepGrid], dtype=np.float32)
    sizes = np.array(n, dtype=np.int32)
    potential_establish = mod.get_function("potential_establish")
    gd = int((n[0]-2)*(n[1]-2)*(n[2]-2)/300)
    if (gd != (n[0]-2)*(n[1]-2)*(n[2]-2)/300):
        gd += 1
    bd = int((n[0]-2)*(n[1]-2)*(n[2]-2)/gd)
    if (bd != (n[0]-2)*(n[1]-2)*(n[2]-2)/gd):
        bd += 1
    print('gd={}, bd = {}'.format(gd, bd))

    prev_phi_gpu = drv.mem_alloc((prev_phi.size) * prev_phi.dtype.itemsize)
    next_phi_gpu = drv.mem_alloc((next_phi.size) * next_phi.dtype.itemsize)
    drv.memcpy_htod(prev_phi_gpu, prev_phi)
    drv.memcpy_htod(next_phi_gpu, next_phi)

    subSum = np.zeros((gd, ), dtype=np.float32)
    ssum = 3 * epsilon
    iter_count = 0
    while ssum > epsilon:
        potential_establish(
            # drv.InOut(prev_phi), drv.InOut(next_phi), drv.InOut(subSum), 
            prev_phi_gpu, next_phi_gpu, drv.InOut(subSum), 
            drv.In(step), drv.In(sizes), 
            block=(bd,1, 1), grid=(gd,1))
        for s in range(subSum.shape[0]):
            subSum[s] = 0.0
        potential_establish(
            # drv.InOut(prev_phi), drv.InOut(next_phi), drv.InOut(subSum), 
            next_phi_gpu, prev_phi_gpu, drv.InOut(subSum), 
            drv.In(step), drv.In(sizes), 
            block=(bd,1, 1), grid=(gd,1))
        ssum = 0
        for s in subSum:
            ssum += s
        if file_name == 'bm_kernel_1.cu':
            prev_sum, next_sum = 0, 0
            for i in range(1, n[0] - 1):
                for j in range(1, n[1] - 1):
                    for k in range(1, n[2] - 1):
                        prev_sum += prev_phi[i][j][k]
                        next_sum += next_phi[i][j][k]
            ssum = np.abs(prev_sum-next_sum)

        iter_count += 2
        # print('{}>{}'.format(ssum, epsilon))
    print('Iteration count = {}'.format(iter_count))
    drv.memcpy_dtoh(prev_phi, prev_phi_gpu)
    drv.memcpy_dtoh(next_phi, next_phi_gpu)

    return (prev_phi, next_phi)


def main(file_name):
    # Границы сетки
    xRange = np.linspace(0, xDimensionGrid, xNumberStepGrid + 1)
    yRange = np.linspace(0, yDimensionGrid, yNumberStepGrid + 1)
    zRange = np.linspace(0, zDimensionGrid, zNumberStepGrid + 1)
    # Граничные условия для потенциала
    n = (xRange.shape[0]+2, yRange.shape[0]+2, zRange.shape[0]+2)
    print('Bounds: {}'.format((xRange.shape[0], yRange.shape[0], zRange.shape[0])))
    prev_phi, next_phi = make_boundary_conditions_for_potentials_2(50, n)
    # for i in range(prev_phi.shape[0]):
    #     print(prev_phi[i])

    n = prev_phi.shape
    ro = np.zeros([n[0], n[1], n[2]], dtype=np.float32)
    start = time.time()
    if (file_name == None):
        potential_establish_method(prev_phi, next_phi, ro, epsilon=1)
    else:
        prev_phi, next_phi = potential_establish_method_cuda(file_name, prev_phi, next_phi, ro, epsilon=1)
    end = time.time()
    print('Elapsed time = {}'.format(end-start))

    # for i in range(next_phi.shape[2]):
    #     print(next_phi[i, 1:-1, 1:-1])
    #     cs=plt.imshow(next_phi[i, 1:-1, 1:-1], interpolation='nearest')
    #     plt.colorbar(cs)
    #     plt.show()

from optparse import OptionParser
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-m", "--method", 
                    dest="method", help="number of method")

    (options, args) = parser.parse_args()

    files = {'0': None, 
            '1': 'bm_kernel_1.cu', 
            '2': 'bm_kernel_2.cu', 
            '3': 'bm_kernel_3.cu', 
            '4': 'kernel_code_one_start_shared_sum.cu', 
            }
    print('\nMethod: {}'.format(files[options.method]))
    # if (files[options.method] == None):
    #     pass
    # else:
    main(files[options.method])


