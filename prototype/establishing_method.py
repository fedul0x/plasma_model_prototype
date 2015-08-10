# -*- coding: utf-8 -*-

import numpy as np
# import pycuda.driver as drv
# import pycuda.autoinit
# from pycuda.compiler import SourceModule

from constant import *


def dif(x, y, z, step):
    return (x - 2 * y + z) / step ** 2


def make_boundary_conditions(phi, n, electron_charge, carbon_charge, helium_charge):
    prev_phi = np.zeros([n[0], n[1], n[2]], dtype=np.float32)
    next_phi = np.zeros([n[0], n[1], n[2]], dtype=np.float32)
    ro = np.zeros([n[0], n[1], n[2]], dtype=np.float32)
    intr = np.linspace(0, phi, n[0])
    for i in range(n[0]):
        for j in range(n[1]):
            for k in range(n[2]):
                if (i * j * k == 0) or (i == n[0] - 1) or (j == n[1] - 1) or (k == n[2] - 1):
                    prev_phi[i][j][k] = intr[-i-1]
                    next_phi[i][j][k] = intr[-i-1]
                    continue
                ro[i][j][k] = \
                    electron_charge[i-1][j-1][k-1] + carbon_charge[i-1][j-1][k-1] + helium_charge[i-1][j-1][k-1]
                prev_phi[i][j][k] = 4*np.pi*ro[i][j][k]*FAKE_TIME_STEP
                next_phi[i][j][k] = 4*np.pi*ro[i][j][k]*FAKE_TIME_STEP
                if i == n[0]-1:
                    prev_phi[i][j][k] = 0
                    next_phi[i][j][k] = 0
                if i == 0:
                    prev_phi[i][j][k] = phi
                    next_phi[i][j][k] = phi

    return (prev_phi, next_phi, ro)


def potential_establish_method(prev_phi, next_phi, ro, epsilon=0.01):
    n = prev_phi.shape
    prev_sum, next_sum = epsilon, 3 * epsilon
    while np.abs(prev_sum - next_sum) > epsilon:
        prev_sum, next_sum = 0, 0
        for i in range(1, n[0] - 1):
            for j in range(1, n[1] - 1):
                for k in range(1, n[2] - 1):
                    next_phi[i][j][k] = FAKE_TIME_STEP * (
                        dif(prev_phi[i + 1][j][k], prev_phi[i][j][k], prev_phi[i - 1][j][k], X_STEP) +
                        dif(prev_phi[i][j + 1][k], prev_phi[i][j][k], prev_phi[i][j - 1][k], Y_STEP) +
                        dif(prev_phi[i][j][k + 1], prev_phi[i][j][k], prev_phi[i][j][k - 1], Z_STEP) +
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
    step = np.array([X_STEP, Y_STEP, Z_STEP], dtype=np.float32)
    sizes = np.array(n, dtype=np.int32)
    potential_establish = mod.get_function("potential_establish")
    gd = int((n[0]-2)*(n[1]-2)*(n[2]-2)/256)
    if (gd != (n[0]-2)*(n[1]-2)*(n[2]-2)/256):
        gd += 1
    bd = int((n[0]-2)*(n[1]-2)*(n[2]-2)/gd)
    if (bd != (n[0]-2)*(n[1]-2)*(n[2]-2)/gd):
        bd += 1
    # print('gd={}, bd = {} ===================='.format(gd, bd))
    # print('sizes={}'.format(sizes))

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
            drv.In(step), drv.In(sizes), np.float32(FAKE_TIME_STEP),
            block=(bd, 1, 1), grid=(gd, 1))
        for s in range(subSum.shape[0]):
            subSum[s] = 0.0
        potential_establish(
            next_phi_gpu, prev_phi_gpu,
            drv.In(ro), drv.InOut(subSum),
            drv.In(step), drv.In(sizes), np.float32(FAKE_TIME_STEP),
            block=(bd, 1, 1), grid=(gd, 1))
        ssum = 0
        for s in subSum:
            ssum += s
    drv.memcpy_dtoh(prev_phi, prev_phi_gpu)
    drv.memcpy_dtoh(next_phi, next_phi_gpu)
    prev_phi_gpu.free()
    next_phi_gpu.free()

    return (prev_phi, next_phi)
