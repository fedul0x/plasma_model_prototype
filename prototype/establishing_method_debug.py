# -*- coding: utf-8 -*-

import os
import time
import shutil
import numpy as np
import matplotlib.pyplot as plt
import pycuda.driver as drv
import pycuda.autoinit
from pycuda.compiler import SourceModule

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
    gd = int((n[0]-2)*(n[1]-2)*(n[2]-2)/256)
    if (gd != (n[0]-2)*(n[1]-2)*(n[2]-2)/256):
        gd += 1
    prev_sum, next_sum = epsilon, 3 * epsilon
    count = 0
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
        prev_phi, next_phi = next_phi, prev_phi
        subSum = np.zeros((gd, ), dtype=np.float32)
        for i in range(1, n[0]-1):
            for j in range(1, n[1]-1):
                for k in range(1, n[2]-1):
                    # prev_phi[i][j][k] = next_phi[i][j][k]
                    # print('{} {} --- {}, {}, {}'.
                        # format(gd, (i*(n[1]-2)*(n[2]-2)+j*(n[2]-2)+k)//256, i, j, k))
                    _i, _j, _k = i-1, j-1, k-1
                    subSum[(_i*(n[1]-2)*(n[2]-2)+_j*(n[2]-2)+_k)//256] += \
                        np.abs(prev_phi[i][j][k] - next_phi[i][j][k])
                    
        # print('{}>{}'.format(np.abs(prev_sum - next_sum), epsilon))
        count += 1
        if count % 10 == 0:
            yield(prev_phi, next_phi, subSum)
    return (prev_phi, next_phi)


def potential_establish_method_cuda(prev_phi, next_phi, ro, epsilon=0.01, delta_t=1.75E-10):
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
    count = 0
    while ssum > epsilon:
        potential_establish(
            prev_phi_gpu, next_phi_gpu,
            drv.In(ro), drv.InOut(subSum),
            drv.In(step), drv.In(sizes), np.float32(FAKE_TIME_STEP),
            block=(bd, 1, 1), grid=(gd, 1))
        for s in range(subSum.shape[0]):
            subSum[s] = 0.0
        count += 1
        potential_establish(
            next_phi_gpu, prev_phi_gpu,
            drv.In(ro), drv.InOut(subSum),
            drv.In(step), drv.In(sizes), np.float32(FAKE_TIME_STEP),
            block=(bd, 1, 1), grid=(gd, 1))
        count += 1
        ssum = 0
        for s in subSum:
            ssum += s
        if count % 10 == 0:
            drv.memcpy_dtoh(prev_phi, prev_phi_gpu)
            drv.memcpy_dtoh(next_phi, next_phi_gpu)
            yield(prev_phi, next_phi, subSum)
    prev_phi_gpu.free()
    next_phi_gpu.free()

    return (prev_phi, next_phi)


def make_dir(prefix, dir_name):
    path = os.path.join(prefix, dir_name)
    try:
        os.mkdir(path)
    except:
        pass
    return path


def make_dir_prefix():
    prefix = './picts/phi_cuda_{}'.format(time.strftime('%Y%m%d%H%M%S'))
    try:
        os.mkdir(prefix)
    except:
        pass
    shutil.copy('./constant.py', prefix)
    return prefix


def make_potential_plot_file(prefix, next_phi_1, next_phi_2, sums_1, sums_2, time):
    directory = make_dir(prefix, 'phi')
    fig = plt.figure()
    ax_1 = fig.add_subplot(221)
    ax_2 = fig.add_subplot(222)
    ax_3 = fig.add_subplot(223)
    ax_4 = fig.add_subplot(224)
    n = next_phi_1.shape
    plot_phi = np.zeros((n[0], n[1]))
    for i in range(n[0]):
        for j in range(n[1]):
            plot_phi[i][j] = next_phi_1[i][j][n[2]//2]
    m1 = ax_1.contourf(plot_phi.T)
    fig.colorbar(m1)
    # ax_1.colorbar()
    for i in range(n[0]):
        for j in range(n[1]):
            plot_phi[i][j] = next_phi_2[i][j][n[2]//2]
    m2 = ax_2.contourf(plot_phi.T)
    fig.colorbar(m2)
    # ax_2.colorbar()
    ax_3.plot(sums_1)
    ax_3.set_title(sum(sums_1))
    ax_4.plot(sums_2)
    ax_4.set_title(sum(sums_2))
    plt.savefig("{}/phi_time={:04d}".format(directory, time))
    plt.clf()
    plt.close()


def main():
    x_range = np.linspace(0, X_DIMENSION_GRID, X_STEP_NUMBER_GRID + 1)
    y_range = np.linspace(0, Y_DIMENSION_GRID, Y_STEP_NUMBER_GRID + 1)
    z_range = np.linspace(0, Z_DIMENSION_GRID, Z_STEP_NUMBER_GRID + 1)
    electron_charge_grid = np.zeros([x_range.shape[0], y_range.shape[0], z_range.shape[0]])
    carbon_charge_grid = np.zeros([x_range.shape[0], y_range.shape[0], z_range.shape[0]])
    helium_charge_grid = np.zeros([x_range.shape[0], y_range.shape[0], z_range.shape[0]])
    for i in range(electron_charge_grid.shape[0]-1):
        for j in range(electron_charge_grid.shape[0]-1):
            for k in range(electron_charge_grid.shape[0]-1):
                electron_charge_grid[i][j][k] = 10E-7
                carbon_charge_grid[i][j][k] = 10E-7
                helium_charge_grid[i][j][k] = 10E-7

    n = (x_range.shape[0]+2, y_range.shape[0]+2, z_range.shape[0]+2)
    phi = POTENTIAL_BOUND_VALUE
    ecg, ccg, hcg = \
        electron_charge_grid, carbon_charge_grid, helium_charge_grid
    prev_phi, next_phi, ro = \
        make_boundary_conditions(phi, n, ecg, ccg, hcg)
    prev_phi_2, next_phi_2 = \
        prev_phi.copy(), next_phi.copy()
        # make_boundary_conditions(phi/2, n, ecg, ccg, hcg)

    generator_1 = potential_establish_method(prev_phi, next_phi, ro, epsilon=ESTABLISHING_METHOD_ACCURACY)
    generator_2 = potential_establish_method_cuda(prev_phi_2, next_phi_2, ro, epsilon=ESTABLISHING_METHOD_ACCURACY)
    # Метод установления
    prefix = make_dir_prefix()
    for i in range(300):
        # next_phi_1, _, sums_1 = next(generator_1)
        # next_phi_2, sums_2 = next_phi_1, sums_1
        next_phi_2, _, sums_2 = next(generator_2)
        next_phi_1, sums_1 = next_phi_2, sums_2
        # print(sums)
        # print(prev_phi, next_phi, sums)
        make_potential_plot_file(prefix, next_phi_1, next_phi_2, sums_1, sums_2, i)
        # prev_phi_2, next_phi_2, sums = next(generator_2)
        # make_potential_plot_file(prefix, next_phi_2, i)


if __name__ == '__main__':
    main()