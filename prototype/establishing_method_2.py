# -*- coding: utf-8 -*-

import numpy as np
import pycuda.driver as drv
import pycuda.autoinit
from pycuda.compiler import SourceModule

from constant import *

def make_boundary_conditions_for_potentials(phi, n, boundary_type='50_50_BOUNDARY'):
    prev_phi = np.zeros([n[0], n[1], n[2]])
    next_phi = np.zeros([n[0], n[1], n[2]])
    middle_phiphi = phi - 0.0000005
    if boundary_type == '50_50_BOUNDARY':
        for i in range(n[0]):
            for j in range(n[1]):
                for k in range(n[2]):
                    if (i * j * k == 0) or (i == n[0] - 1) or (j == n[1] - 1) or (k == n[2] - 1):
                        prev_phi[i][j][k] = phi
                        next_phi[i][j][k] = phi
                    else:
                        prev_phi[i][j][k] = middle_phiphi
                        next_phi[i][j][k] = middle_phiphi
    if boundary_type == '50_0_INTERPOLATION_BOUNDARY':
        intr = np.linspace(0, 50, n[0])
        for j in range(n[1]):
            for k in range(n[2]):
                for i in range(n[0]):
                        prev_phi[i][j][k] = intr[-i-1] # + np.random.random()*2
                        next_phi[i][j][k] = intr[-i-1] # + np.random.random()*2
    return (prev_phi, next_phi)


def make_boundary_conditions_for_potentials_2(phi, n, chargeGridElectron, chargeGridCarbon, chargeGridHelium, boundary_type='50_50_BOUNDARY'):
    prev_phi = np.zeros([n[0], n[1], n[2]], dtype=np.float32)
    next_phi = np.zeros([n[0], n[1], n[2]], dtype=np.float32)
    middle_phiphi = phi - 0.0000005
    n2 = chargeGridElectron.shape
    if boundary_type == '50_50_BOUNDARY':
        for i in range(n[0]):
            for j in range(n[1]):
                for k in range(n[2]):
                    if i-1>=n2[0] or j-1>=n2[1] or k-1>=n2[2] or i-1<0 or j-1<0 or k-1<0:
                        prev_phi[i][j][k] = phi
                        # next_phi[i][j][k] = 0
                        next_phi[i][j][k] = phi
                    else:
                        ro = chargeGridElectron[i-1][j-1][k-1] + chargeGridCarbon[i-1][j-1][k-1] + chargeGridHelium[i-1][j-1][k-1]
                        prev_phi[i][j][k] = 4*np.pi*ro*deltaT_
                        next_phi[i][j][k] = 0
    if boundary_type == '50_0_BOUNDARY':
        ro = chargeGridElectron[0][0][0] + chargeGridCarbon[0][0][0] + chargeGridHelium[0][0][0]
        for i in range(n[0]):
            for j in range(n[1]):
                for k in range(n[2]):
                    if not (i-1>=n2[0]-1 or j-1>=n2[1]-1 or k-1>=n2[2]-1 or i-1<0 or j-1<0 or k-1<0):
                        # pass
                        ro = chargeGridElectron[i][j][k] + chargeGridCarbon[i][j][k] + chargeGridHelium[i][j][k]
                    prev_phi[i][j][k] = 4*np.pi*ro*deltaT_
                    next_phi[i][j][k] = 4*np.pi*ro*deltaT_
                    if i == n[0]-1: 
                        prev_phi[i][j][k] = 0
                        next_phi[i][j][k] = 0
                    if i == 0: 
                        prev_phi[i][j][k] = phi
                        next_phi[i][j][k] = phi

    # print('THISSSS')
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


def potential_establish_method_cuda(prev_phi, next_phi, ro, epsilon=0.01):

    f = open('./kernel_1.cu', 'r')
    mod = SourceModule("".join(f.readlines()))
    n = prev_phi.shape
    step = np.array([xStepGrid, yStepGrid, zStepGrid], dtype=np.float32)
    sizes = np.array(n, dtype=np.int32)
    potential_establish = mod.get_function("potential_establish")
    gd = int((n[0]-2)*(n[1]-2)*(n[2]-2)/256)
    if (gd != (n[0]-2)*(n[1]-2)*(n[2]-2)/256):
        gd += 1
    bd = int((n[0]-2)*(n[1]-2)*(n[2]-2)/gd)
    if (bd != (n[0]-2)*(n[1]-2)*(n[2]-2)/gd):
        bd += 1
    print('gd={}, bd = {} ===================='.format(gd, bd))
    print('sizes={}'.format(sizes))

    prev_phi_gpu = drv.mem_alloc((prev_phi.size) * prev_phi.dtype.itemsize)
    next_phi_gpu = drv.mem_alloc((next_phi.size) * next_phi.dtype.itemsize)
    drv.memcpy_htod(prev_phi_gpu, prev_phi)
    drv.memcpy_htod(next_phi_gpu, next_phi)

    subSum = np.zeros((gd, ), dtype=np.float32)
    ssum = 4 * epsilon
    while ssum > epsilon:
        potential_establish(
            prev_phi_gpu, next_phi_gpu, 
            drv.In(ro), drv.InOut(subSum), 
            drv.In(step), drv.In(sizes), 
            block=(bd,1, 1), grid=(gd,1))
        for s in range(subSum.shape[0]):
            subSum[s] = 0.0
        potential_establish(
            next_phi_gpu, prev_phi_gpu, 
            drv.In(ro), drv.InOut(subSum), 
            drv.In(step), drv.In(sizes), 
            block=(bd,1, 1), grid=(gd,1))
        ssum = 0
        for s in subSum:
            ssum += s
        # print('{}>{}'.format(ssum, epsilon))

    # for i in range(flags.shape[0]):
    #     # cs = plt.contourf(flags[i])
    #     cs=plt.imshow(flags[i], interpolation='nearest')
    #     plt.colorbar(cs)
    #     plt.show()
    drv.memcpy_dtoh(prev_phi, prev_phi_gpu)
    drv.memcpy_dtoh(next_phi, next_phi_gpu)
    prev_phi_gpu.free()
    next_phi_gpu.free()

    return (prev_phi, next_phi)
