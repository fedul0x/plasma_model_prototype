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
                            prev_phi[i][j][k] = 4*np.pi*ro*deltaT_
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

    mod = SourceModule("""

    __device__ __host__ float dif(float x, float y, float z, float step){
        return (x - 2 * y + z) / step / step;
    }

    //TODO change step and matrix_size to vector
    __global__ void potential_establish(float* prev_phi, float* next_phi, float* step, int* matrix_size )
    {
        const int x_dim = matrix_size[0];
        const int y_dim = matrix_size[1];
        const int z_dim = matrix_size[2];

        const int index = blockIdx.x*blockDim.x + threadIdx.x;
        if (index >=(x_dim-2)*(y_dim-2)*(z_dim-2)) 
            return;
        const int i = (int) (index/(z_dim-2)/(y_dim-2))+1;
        const int j = ((int) (index/(z_dim-2)))%(y_dim-2)+1;
        const int k = index%(z_dim-2)+1;
        // const int i = (int) (threadIdx.x/(z_dim-2)/(y_dim-2))+1;
        // const int j = ((int) (threadIdx.x/(z_dim-2)))%(y_dim-2)+1;
        // const int k = threadIdx.x%(z_dim-2)+1;

        const float x_step = step[0];
        const float y_step = step[1];
        const float z_step = step[2];

        const int ijk = i*z_dim*y_dim + j*z_dim + k;
        const int im1jk = (i-1)*z_dim*y_dim + j*z_dim + k;
        const int ijm1k = i*z_dim*y_dim + (j-1)*z_dim + k;
        const int ijkm1 = i*z_dim*y_dim + j*z_dim + k - 1;
        //const int ijkm1 = ijk - 1;
        const int ip1jk = (i+1)*z_dim*y_dim + j*z_dim + k;
        const int ijp1k = i*z_dim*y_dim + (j+1)*z_dim + k;
        const int ijkp1 = i*z_dim*y_dim + j*z_dim + k + 1;
        //const int ijkp1 = ijk + 1;

        //for(int l=0; l <= 1000000; l++) 
        float tmp = prev_phi[ijk];
        float ro = 1E-17;
        //for(int l=0; l <= 500; l++) 
        {
            prev_phi[ijk] = next_phi[ijk];
            //TODO rewrite PI and deltaT_
            //next_phi[i][j][k] = 1.75E-10 * (dif(prev_phi[i + 1][j][k], prev_phi[i][j][k], prev_phi[i - 1][j][k], x_step) + dif(prev_phi[i][j + 1][k], prev_phi[i][j][k], prev_phi[i][j - 1][k], y_step) + dif(prev_phi[i][j][k + 1], prev_phi[i][j][k], prev_phi[i][j][k - 1], z_step) + 4 * 3.1415926 * ro[i][j][k]) + prev_phi[i][j][k];
            float tmp = 1.75E-10 * (
            dif(prev_phi[ip1jk], prev_phi[ijk], prev_phi[im1jk], x_step) +
            dif(prev_phi[ijp1k], prev_phi[ijk], prev_phi[ijm1k], y_step) +
            dif(prev_phi[ijkp1], prev_phi[ijk], prev_phi[ijkm1], z_step) +
            4 * 3.1415926 * ro) + prev_phi[ijk];
            //tmp = dif(prev_phi[ip1jk], prev_phi[ijk], prev_phi[im1jk], x_step);
            __syncthreads();
            //TODO change this append to another where not depend to next 
            //next_phi[ijk] = tmp;
            //next_phi[ijk] = ijk + 1000*ijkm1 + 1000000*ijkp1;
            //next_phi[ijk] = 1000*ijk + 1000000*ijkm1 + ijkp1;
            next_phi[ijk] = 1000*next_phi[ijk] + 1000000*next_phi[ijkm1]
             + next_phi[ijkp1];
            //next_phi[ijk] = 100*ijk + 10000*ijm1k + ijp1k;
            //next_phi[ijk] = next_phi[ijk] + 1;
        }
    } 
    """)
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
    print(gd, bd)
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
        for kk in range(prev_phi.shape[0]):
            print(next_phi[kk])
            # plt.contourf(next_phi[kk])
            # plt.show()
        print('{}>{}'.format(np.abs(prev_sum - next_sum), epsilon))
        # break
        # print('===================')
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
    n = prev_phi.shape
    ro = np.zeros([n[0], n[1], n[2]], dtype=np.float32)
    # for i in range(1, n[0] - 1):
    #     for j in range(1, n[1] - 1):
    #         for k in range(1, n[2] - 1):
    #             ro[i][j][k] = chargeGridElectron[i-1][j-1][k-1] + chargeGridCarbon[i-1][j-1][k-1] + chargeGridHelium[i-1][j-1][k-1]
    print('THIS IS START OF potential_establish_method_2')
    prev_phi, next_phi = potential_establish_method_2(prev_phi, next_phi, ro, epsilon=0.01)
    print('THIS IS END OF potential_establish_method_2')

    for i in range(next_phi.shape[0]):
        print(next_phi[i])
        plt.contourf(next_phi[i])
        plt.show()

if __name__ == '__main__':
    main()
