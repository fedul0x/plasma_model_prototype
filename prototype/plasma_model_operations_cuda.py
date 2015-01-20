# -*- coding: utf-8 -*-
import scipy
import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from constant import *
from plot import *


# import pdb

''' pseuso time '''
def p_time(time):
    return time % (in_memory_time)

def p_next_time(time):
    return (time + 1) % (in_memory_time)

def p_prev_time(time):
    return (time - 1) % (in_memory_time)

''' Maxwell distribution randomizer '''
class MWranomizer(scipy.stats.rv_continuous):
    def __init__(self):
        scipy.stats.rv_continuous.__init__(self)
        self.a = 0

    def _pdf(self, v, n, m, T):
        Boltzmann_constant = 1.380648813E-23
        k = Boltzmann_constant
        return n*(m/(2*np.pi*k*T))**(3/2)*np.exp(-m*v*v/(2*k*T))*4*np.pi*v**2
        # return n*np.exp(v)*4*np.pi*v**2

def get_component(particle, b=0, n=3):
    return [particle[i] for i in range(b, b+n)]

def getSpeedProjection(particle):
    # speed = particle[5]
    # devider = np.sqrt(3)
    # return (speed*Kvu￼1/devider, speed*Kvu2/devider, speed*Kvu2/devider)
    return (particle[5], particle[6], particle[7])

def spreadParticle(center, steps, speedCenter, number):
    x, y, z, = center
    xStepGrid, yStepGrid, zStepGrid = steps
    xBig, yBig, zBig, speed = 0, 0, 0, 0
    maxSpeed, deltaSpeed = speedCenter
    for _ in range(number):
        xCurrent = x + rand() * xStepGrid
        yCurrent = y + rand() * yStepGrid
        zCurrent = z + rand() * zStepGrid
        xBig += xCurrent
        yBig += yCurrent
        zBig += zCurrent
        speed += (maxSpeed + deltaSpeed) - 2*rand()*deltaSpeed
    # return xBig/number/Ml, yBig/number/Ml, zBig/number/Ml, speed/number
    return xBig/number, yBig/number, zBig/number, speed/number

def spreadCharge(center, steps, bigCenter, charge):
    x, y, z, = center
    i, j, k = 0, 0, 0
    xStepGrid, yStepGrid, zStepGrid = steps
    xBig, yBig, zBig = bigCenter
    v = xStepGrid * yStepGrid * zStepGrid
    neighborX, neighborY, neighborZ = xBig - x, yBig - y, zBig - z
    oppositeX, oppositeY, oppositeZ = xStepGrid - xBig + x, yStepGrid - yBig + y, zStepGrid - zBig + z

    res = []
    res += [(i, j, k, charge * oppositeX * oppositeY * oppositeZ / v)]
    res += [(i + 1, j, k, charge * neighborX * oppositeY * oppositeZ / v)]
    res += [(i, j + 1, k, charge * oppositeX * neighborY * oppositeZ / v)]
    res += [(i, j, k + 1, charge * oppositeX * oppositeY * neighborZ / v)]
    res += [(i + 1, j + 1, k, charge * neighborX * neighborY * oppositeZ / v)]
    res += [(i, j + 1, k + 1, charge * oppositeX * neighborY * neighborZ / v)]
    res += [(i + 1, j, k + 1, charge * neighborX * oppositeY * neighborZ / v)]
    res += [(i + 1, j + 1, k + 1, charge * neighborX * neighborY * neighborZ / v)]

    return res

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
        for i in range(n[0]):
            for j in range(n[1]):
                for k in range(n[2]):
                    ro = chargeGridElectron[i-1][j-1][k-1] + chargeGridCarbon[i-1][j-1][k-1] + chargeGridHelium[i-1][j-1][k-1]
                    prev_phi[i][j][k] = 4*np.pi*ro*deltaT_
                    next_phi[i][j][k] = 4*np.pi*ro*deltaT_
                    if i == n2[0]-1: 
                        prev_phi[i][j][k] = 0
                        next_phi[i][j][k] = 0
                    if i == 0: 
                        prev_phi[i][j][k] = phi
                        next_phi[i][j][k] = phi

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
    import pycuda
    import pycuda.driver as drv
    import pycuda.autoinit
    import numpy
    from pycuda.compiler import SourceModule

    mod = SourceModule("""

    __device__ __host__ float dif(float x, float y, float z, float step){
        return (x - 2 * y + z) / step / step;
    }

    //__global__ void potential_establish(float*** prev_phi, float*** next_phi, float*** ro, float* step, float epsilon)
    //TODO change step and matrix_size to vector
    __global__ void potential_establish(float* prev_phi, float* next_phi, float* ro, float* step, int* matrix_size )
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
        //for(int l=0; l <= 500; l++) 
        {
            prev_phi[ijk] = next_phi[ijk];
            //TODO rewrite PI and deltaT_
            //next_phi[i][j][k] = 1.75E-10 * (dif(prev_phi[i + 1][j][k], prev_phi[i][j][k], prev_phi[i - 1][j][k], x_step) + dif(prev_phi[i][j + 1][k], prev_phi[i][j][k], prev_phi[i][j - 1][k], y_step) + dif(prev_phi[i][j][k + 1], prev_phi[i][j][k], prev_phi[i][j][k - 1], z_step) + 4 * 3.1415926 * ro[i][j][k]) + prev_phi[i][j][k];
            float tmp = 1.75E-10 * (
            dif(prev_phi[ip1jk], prev_phi[ijk], prev_phi[im1jk], x_step) +
            dif(prev_phi[ijp1k], prev_phi[ijk], prev_phi[ijm1k], y_step) +
            dif(prev_phi[ijkp1], prev_phi[ijk], prev_phi[ijkm1], z_step) +
            4 * 3.1415926 * ro[ijk]) + prev_phi[ijk];
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
    # potential_establish(
    #     drv.Out(prev_phi), drv.Out(next_phi), drv.In(ro), drv.In(step), drv.In(epsilon),
    #     block=(n[1]-2,n[2]-2,1), grid=(n[0]-2,1))
    # print(next_phi)
    # print('===================')
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
            drv.InOut(prev_phi), drv.InOut(next_phi), drv.In(ro), drv.In(step), drv.In(sizes), 
            # block=((n[0]-2)*(n[1]-2)*(n[2]-2),1, 1), grid=(1,1))
            block=(bd,1, 1), grid=(gd,1))
            # block=(512,1, 1), grid=(1,1))
            # block=(n[0],n[1],n[2]), grid=(1,1))
            # block=(n[1]-2,n[2]-2,1), grid=(n[0]-2,1))
            # block=(400,1,1), grid=(1,1))
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


def spreadTension(center, steps, bigCenter, intensity):
    x, y, z, = center
    xStepGrid, yStepGrid, zStepGrid = steps
    xBig, yBig, zBig = bigCenter
    i, j, k = int(xBig/xStepGrid), int(yBig/yStepGrid), int(zBig/zStepGrid)
    v = xStepGrid * yStepGrid * zStepGrid
    neighborX, neighborY, neighborZ = xBig - x, yBig - y, zBig - z
    oppositeX, oppositeY, oppositeZ = xStepGrid - xBig + x, yStepGrid - yBig + y, zStepGrid - zBig + z
    xTension = intensity[i][j][k][0] * (xStepGrid-xBig+x) * (yStepGrid-yBig+y) * (zStepGrid-zBig+z) / v +\
        intensity[i+1][j][k][0] * (neighborX)*(yStepGrid-yBig+y)*(zStepGrid-zBig+z) / v +\
        intensity[i][j+1][k][0] * (xStepGrid-xBig+x)*(neighborY)*(zStepGrid-zBig+z) / v +\
        intensity[i][j][k+1][0] * (xStepGrid-xBig+x)*(yStepGrid-yBig+y)*neighborZ / v +\
        intensity[i+1][j+1][k][0] * (neighborX)*(neighborY)*(zStepGrid-zBig+z) / v +\
        intensity[i][j+1][k+1][0] * (xStepGrid-xBig+x)*(neighborY)*neighborZ / v +\
        intensity[i+1][j][k+1][0] * (neighborX)*(yStepGrid-yBig+y)*neighborZ / v +\
        intensity[i+1][j+1][k+1][0] * (neighborX)*(neighborY)*neighborZ / v

    yTension = intensity[i][j][k][1] * (xStepGrid-xBig+x) * (yStepGrid-yBig+y) * (zStepGrid-zBig+z) / v +\
        intensity[i+1][j][k][1] * (neighborX)*(yStepGrid-yBig+y)*(zStepGrid-zBig+z) / (xStepGrid*yStepGrid*zStepGrid) +\
        intensity[i][j+1][k][1] * (xStepGrid-xBig+x)*(neighborY)*(zStepGrid-zBig+z) / (xStepGrid*yStepGrid*zStepGrid) +\
        intensity[i][j][k+1][1] * (xStepGrid-xBig+x)*(yStepGrid-yBig+y)*neighborZ / (xStepGrid*yStepGrid*zStepGrid) +\
        intensity[i+1][j+1][k][1] * (neighborX)*(neighborY)*(zStepGrid-zBig+z) / (xStepGrid*yStepGrid*zStepGrid) +\
        intensity[i][j+1][k+1][1] * (xStepGrid-xBig+x)*(neighborY)*neighborZ / (xStepGrid*yStepGrid*zStepGrid) +\
        intensity[i+1][j][k+1][1] * (neighborX)*(yStepGrid-yBig+y)*neighborZ / (xStepGrid*yStepGrid*zStepGrid) +\
        intensity[i+1][j+1][k+1][1] * (neighborX)*(neighborY)*neighborZ / (xStepGrid*yStepGrid*zStepGrid)

    zTension = intensity[i][j][k][2] * (xStepGrid-xBig+x) * (yStepGrid-yBig+y) * (zStepGrid-zBig+z) / v +\
        intensity[i+1][j][k][2] * (neighborX)*(yStepGrid-yBig+y)*(zStepGrid-zBig+z) / (xStepGrid*yStepGrid*zStepGrid) +\
        intensity[i][j+1][k][2] * (xStepGrid-xBig+x)*(neighborY)*(zStepGrid-zBig+z) / (xStepGrid*yStepGrid*zStepGrid) +\
        intensity[i][j][k+1][2] * (xStepGrid-xBig+x)*(yStepGrid-yBig+y)*neighborZ / (xStepGrid*yStepGrid*zStepGrid) +\
        intensity[i+1][j+1][k][2] * (neighborX)*(neighborY)*(zStepGrid-zBig+z) / (xStepGrid*yStepGrid*zStepGrid) +\
        intensity[i][j+1][k+1][2] * (xStepGrid-xBig+x)*(neighborY)*neighborZ / (xStepGrid*yStepGrid*zStepGrid) +\
        intensity[i+1][j][k+1][2] * (neighborX)*(yStepGrid-yBig+y)*neighborZ / (xStepGrid*yStepGrid*zStepGrid) +\
        intensity[i+1][j+1][k+1][2] * (neighborX)*(neighborY)*neighborZ / (xStepGrid*yStepGrid*zStepGrid)
    return (xTension, yTension, zTension)

def dif(x, y, z, step):
    return (x - 2 * y + z) / step ** 2




def main(prefix):
    # Границы сетки
    xRange = np.linspace(0, xDimensionGrid, xNumberStepGrid + 1)
    yRange = np.linspace(0, yDimensionGrid, yNumberStepGrid + 1)
    zRange = np.linspace(0, zDimensionGrid, zNumberStepGrid + 1)
    n = (xRange.shape[0] - 1) * (yRange.shape[0] - 1) * (zRange.shape[0] - 1) 
    # time, i, x, y, z, radius, charge, speedx, speedy, speedz 
    positionElectron = np.empty([in_memory_time, n, 6+2])
    positionCarbon = np.empty([in_memory_time, int(n/(xRange.shape[0] - 1)), 6+2])
    positionHelium = np.empty([in_memory_time, n, 6+2])
    # Заряды в узлах сетки
    # chargeGridElectron = np.zeros([xRange.shape[0], yRange.shape[0], zRange.shape[0]])
    # chargeGridCarbon = np.zeros([xRange.shape[0], yRange.shape[0], zRange.shape[0]])
    # chargeGridHelium = np.zeros([xRange.shape[0], yRange.shape[0], zRange.shape[0]])
    # Сразу крупные частицы
    num = 0 # номер частицы
    numCarbon = 0 # номер углерода
    time = 0 #позиция по временной шкале с учетом цикличности заполнения
    for x, i in zip(xRange[:-1], range(xRange.shape[0] - 1)):
        for y, j in zip(yRange[:-1], range(yRange.shape[0] - 1)):
            for z, k in zip(zRange[:-1], range(zRange.shape[0] - 1)):
                # Распределение электронов
                xBig, yBig, zBig, speed = spreadParticle((x, y, z), (xStepGrid, yStepGrid, zStepGrid), (maxSpeedElectron, deltaSpeedElectron), numbersElectron)
                speed = speed/Mnue
                # холодный страт электронов
                speed = 0.0
                radiusBigElectron = (hi * radiusElectron ** 3 * numbersElectron) ** (1.0 / 3.0)
                chargeBigElectron = chargeElectron * numbersElectron
                for v, l in zip([xBig, yBig, zBig, radiusBigElectron, chargeBigElectron, speed, speed, speed], range(positionElectron.shape[2])):
                    positionElectron[time][num][l] = v

                # Распределение углерода
                if i == 0:
                    xBig, yBig, zBig, speed = spreadParticle((0, y, z), (xStepGrid, yStepGrid, zStepGrid), (maxSpeedCarbon, deltaSpeedCarbon), numbersCarbon)
                    speed = speed/Mnuc
                    # TODO горячий страт
                    speed = 0.0
                    radiusBigCarbon = (hi * radiusCarbon ** 3 * numbersCarbon) ** (1.0 / 3.0) * 10.0**6
                    chargeBigCarbon = chargeCarbon * numbersCarbon
                    for v, l in zip([0, yBig, zBig, radiusBigCarbon, chargeBigCarbon, speed, speed, speed], range(positionCarbon.shape[2])):
                        positionCarbon[time][numCarbon][l] = v
                    numCarbon += 1

                # Распределение гелия
                xBig, yBig, zBig, speed = spreadParticle((x, y, z), (xStepGrid, yStepGrid, zStepGrid), (maxSpeedHelium, deltaSpeedHelium), numbersHelium)
                speed = speed/Mnuh
                # TODO горячий страт
                speed = 0.0
                radiusBigHelium = (hi * radiusHelium ** 3 * numbersHelium) ** (1.0 / 3.0)
                chargeBigHelium = chargeHelium * numbersHelium
                for v, l in zip([xBig, yBig, zBig, radiusBigHelium, chargeBigHelium, speed, speed, speed], range(positionHelium.shape[2])):
                    positionHelium[time][num][l] = v

                num += 1

    
    # MODELING CYCLE BEGIN
    num = 0 # номер частицы
    time = 0 # абсолютная позиция по временной шкале
    crashesElectron = []
    lastCrashesElectron = []
    crashesHelium = []
    lastCrashesHelium = []
    # Для итоговых графиков
    listen_particles = [0, int(positionCarbon.shape[1]/2), positionCarbon.shape[1]-1]
    plot_data = [[] for _ in listen_particles]
    try:
        while (time < modeling_time):
            curr_time = p_time(time)
            # Заряд в узлах
            chargeGridElectron = np.zeros([xRange.shape[0], yRange.shape[0], zRange.shape[0]])
            chargeGridCarbon = np.zeros([xRange.shape[0], yRange.shape[0], zRange.shape[0]])
            chargeGridHelium = np.zeros([xRange.shape[0], yRange.shape[0], zRange.shape[0]])
            positions = [positionElectron, positionCarbon, positionHelium]
            grids = [chargeGridElectron, chargeGridCarbon, chargeGridHelium]
            for grid, position in zip(grids, positions):
                for num in range(position.shape[1]):
                    xBig, yBig, zBig, _, charge, _, _, _ = position[curr_time][num]
                    # if xBig/xStepGrid == np.NaN or yBig/yStepGrid == np.NaN or zBig/zStepGrid == np.NaN:
                    #     pdb.set_trace()
                    i, j, k = int(xBig/xStepGrid), int(yBig/yStepGrid), int(zBig/zStepGrid)
                    x, y, z = i*xStepGrid, j*yStepGrid, k*zStepGrid
                    patch = spreadCharge((x, y, z), (xStepGrid, yStepGrid, zStepGrid), (xBig, yBig, zBig), charge)
                    for p in patch:
                        # # print(p, i, j, k)
                        # print(xBig, yBig, zBig)
                        # # print(i, j, k)
                        grid[i+p[0]][j+p[1]][k+p[2]] += p[3]
            # Граничные условия для потенциала
            n = (xRange.shape[0]+2, yRange.shape[0]+2, zRange.shape[0]+2)
            # prev_phi, next_phi = make_boundary_conditions_for_potentials(50, n)
            prev_phi, next_phi = make_boundary_conditions_for_potentials_2(50, n, chargeGridElectron, chargeGridCarbon, chargeGridHelium)
            # prev_phi, next_phi = make_boundary_conditions_for_potentials(50, n, boundary_type='50_0_INTERPOLATION_BOUNDARY')
            # Метод установления
            # n = (chargeGridElectron.shape[0]+2, chargeGridElectron.shape[1]+2, chargeGridElectron.shape[2]+2)
            # n = (prev_phi.shape[0], prev_phi.shape[1], prev_phi.shape[2])
            n = prev_phi.shape
            ro = np.zeros([n[0], n[1], n[2]], dtype=np.float32)
            for i in range(1, n[0] - 1):
                for j in range(1, n[1] - 1):
                    for k in range(1, n[2] - 1):
                        ro[i][j][k] = chargeGridElectron[i-1][j-1][k-1] + chargeGridCarbon[i-1][j-1][k-1] + chargeGridHelium[i-1][j-1][k-1]
            # print(prev_phi)
            # print(next_phi)
            print('THIS IS START OF potential_establish_method_2')
            prev_phi, next_phi = potential_establish_method_2(prev_phi, next_phi, ro, epsilon=0.01)
            print('THIS IS END OF potential_establish_method_2')
            # return

            # Расчет напряженности
            n = next_phi.shape
            intensity = np.zeros([n[0]-2, n[1]-2, n[2]-2, 3])
            # pdb.set_trace()
            inten = [[], [], []]
            for i in range(1, n[0]-1):
                for j in range(1, n[1]-1):
                    for k in range(1, n[2]-1):
                        # print(xStepGrid, yStepGrid, zStepGrid)
                        intensity[i-1][j-1][k-1][0] = (next_phi[i - 1][j][k] - next_phi[i + 1][j][k]) / 2 / xStepGrid
                        intensity[i-1][j-1][k-1][1] = (next_phi[i][j - 1][k] - next_phi[i][j + 1][k]) / 2 / yStepGrid
                        intensity[i-1][j-1][k-1][2] = (next_phi[i][j][k - 1] - next_phi[i][j][k + 1]) / 2 / zStepGrid
                        inten[0] += [(intensity[i-1][j-1][k-1][0], next_phi[i - 1][j][k], next_phi[i + 1][j][k])]
                        inten[1] += [(intensity[i-1][j-1][k-1][1], next_phi[i][j - 1][k], next_phi[i][j + 1][k])]
                        inten[2] += [(intensity[i-1][j-1][k-1][2], next_phi[i][j][k - 1], next_phi[i][j][k + 1])]

            # Расчет напряженности действующей на частицу
            tensionCorpusculActingElectron = np.empty([positionElectron.shape[1], 3])
            tensionCorpusculActingCarbon = np.empty([positionCarbon.shape[1], 3])
            tensionCorpusculActingHelium = np.empty([positionHelium.shape[1], 3])
            particles = [positionElectron, positionCarbon, positionHelium]
            tensions = [tensionCorpusculActingElectron, tensionCorpusculActingCarbon, tensionCorpusculActingHelium]
            for position, p in zip(particles, range(len(tensions))):
                n = position.shape
                for num in range(n[1]):
                    xBig, yBig, zBig = get_component(position[curr_time][num])
                    # xBig, yBig, zBig = position[curr_time][num][0], position[curr_time][num][1], position[curr_time][num][2]
                    i, j, k = int(xBig/xStepGrid), int(yBig/yStepGrid), int(zBig/zStepGrid)
                    x, y, z = i*xStepGrid, j*yStepGrid, k*zStepGrid
                    tension = spreadTension((x, y, z), (xStepGrid, yStepGrid, zStepGrid), (xBig, yBig, zBig), intensity)
                    for t in range(3):
                        tensions[p][num][t] = tension[t]


            # # Стокновения углерода и электронов
            # pe = positionElectron
            # pc = positionCarbon
            # lastCrashesElectron = crashesElectron
            # crashesElectron = []

            # for e in range(pe.shape[1]):
            #     for c in range(pc.shape[1]):
            #         # print("{} <= {}".format(np.sqrt((pe[curr_time][e][0] - pc[curr_time][c][0] )**2 + (pe[curr_time][e][1] - pc[curr_time][c][1] )**2 + (pe[curr_time][e][2] - pc[curr_time][c][2] )**2),  pe[curr_time][e][3] + pc[curr_time][c][3]))
            #         # print("{} > {}".format(np.sqrt((pe[curr_time][e][0] - pc[curr_time][c][0] )**2 + (pe[curr_time][e][1] - pc[curr_time][c][1] )**2 + (pe[curr_time][e][2] - pc[curr_time][c][2] )**2), np.abs(pe[curr_time][e][3] - pc[curr_time][c][3])))
            #         if np.sqrt((pe[curr_time][e][0] - pc[curr_time][c][0] )**2 + (pe[curr_time][e][1] - pc[curr_time][c][1] )**2 + (pe[curr_time][e][2] - pc[curr_time][c][2] )**2)  <= pe[curr_time][e][3] + pc[curr_time][c][3]  and \
            #             np.sqrt((pe[curr_time][e][0] - pc[curr_time][c][0] )**2 + (pe[curr_time][e][1] - pc[curr_time][c][1] )**2 + (pe[curr_time][e][2] - pc[curr_time][c][2] )**2) > np.abs(pe[curr_time][e][3] - pc[curr_time][c][3]):

            #             crashesElectron += [pe[curr_time][e].copy()]
            #             # print("++++++++++++ELECTRON+++++++++++++++++++++++")
            #             # print(pe[curr_time][e])
            #             # print("++++++++++++++++++++++++++++++++++++")
            #             spd1 = [(2 * massElectron * pe[curr_time][e][p]  + pc[curr_time][c][p]* (massCarbon - massElectron))/(massCarbon + massElectron) for p in range(5, 8)]
            #             spd2 = [(2 * massCarbon * pc[curr_time][c][p]  + pe[curr_time][e][p]* (massElectron - massCarbon))/(massCarbon + massElectron) for p in range(5, 8)]

            #             pc[curr_time][c][5], pc[curr_time][c][6], pc[curr_time][c][7] = spd1[0], spd1[1], spd1[2]
            #             pe[curr_time][e][5], pe[curr_time][e][6], pe[curr_time][e][7] = spd2[0], spd2[1], spd2[2]
            #             # print(pe[curr_time][e])
            #             # print("+++++++++++++++ELECTRON+END+++++++++++++++++++++")
            #             # print(crashesElectron)
            #             # print("+++++++++++++++CARBON+END+++++++++++++++++++++")
            #             # print(pc[curr_time][c])
            #             # print("+++++++++++++++ELECTRON+END+++++++++++++++++++++")
            # ph = positionHelium
            # pc = positionCarbon
            # lastCrashesHelium = crashesHelium
            # crashesHelium = []
            # for h in range(ph.shape[1]):
            #     for c in range(pc.shape[1]):
            #         # print("{} <= {}".format(np.sqrt((ph[curr_time][h][0] - pc[curr_time][c][0] )**2 + (ph[curr_time][h][1] - pc[curr_time][c][1] )**2 + (ph[curr_time][h][2] - pc[curr_time][c][2] )**2),  ph[curr_time][h][3] + pc[curr_time][c][3]))
            #         # print("{} > {}".format(np.sqrt((ph[curr_time][h][0] - pc[curr_time][c][0] )**2 + (ph[curr_time][h][1] - pc[curr_time][c][1] )**2 + (ph[curr_time][h][2] - pc[curr_time][c][2] )**2), np.abs(ph[curr_time][h][3] - pc[curr_time][c][3])))
            #         if np.sqrt((ph[curr_time][h][0] - pc[curr_time][c][0] )**2 + (ph[curr_time][h][1] - pc[curr_time][c][1] )**2 + (ph[curr_time][h][2] - pc[curr_time][c][2] )**2)  <= ph[curr_time][h][3] + pc[curr_time][c][3] and \
            #             np.sqrt((ph[curr_time][h][0] - pc[curr_time][c][0] )**2 + (ph[curr_time][h][1] - pc[curr_time][c][1] )**2 + (ph[curr_time][h][2] - pc[curr_time][c][2] )**2) > np.abs(ph[curr_time][h][3] - pc[curr_time][c][3]):

            #             crashesHelium += [ph[curr_time][h][:]]
            #             print("+++++++++++++++HELIUM+++++++++++++++++++++")
            #             print(ph[curr_time][h][:])
            #             print("++++++++++++++++++++++++++++++++++++")
            #             spd1 = [(2 * massHelium * ph[curr_time][h][p]  + pc[curr_time][c][p]* (massCarbon - massHelium))/(massCarbon + massHelium) for p in range(5, 8)]
            #             spd2 = [(2 * massCarbon * pc[curr_time][c][p]  + ph[curr_time][h][p]* (massHelium - massCarbon))/(massCarbon + massHelium) for p in range(5, 8)]

            #             pc[curr_time][c][5], pc[curr_time][c][6], pc[curr_time][c][7] = spd1[0], spd1[1], spd1[2]
            #             ph[curr_time][h][5], ph[curr_time][h][6], ph[curr_time][h][7] = spd2[0], spd2[1], spd2[2]
            #             print(ph[curr_time][h][:])
            #             print("+++++++++++++++HELIUM+END+++++++++++++++++++++")
            #             print(crashesHelium)
            #             print("+++++++++++++++HELIUM+END+++++++++++++++++++++")
            print('THIS IS START OF RKF')
            from scipy.integrate import odeint
            def f(y, t):
                # delta = 1.0
                # Z = 1.0
                # epsilonC = 1.0
                r, v, E, delta, Z, epsilonC = y[0], y[1], y[2], y[3], y[4], y[5]
                f = [np.sqrt(delta)*v, np.sqrt(delta)*Z/2/epsilonC*E]
                return f

            t = np.linspace(0, 2*deltaT, 20)
            curr_time = p_time(time)
            for num in range(positionCarbon.shape[1]):
                xBig, yBig, zBig = get_component(positionCarbon[curr_time][num])
                # xBig, yBig, zBig, _, _, _, _, _ = positionCarbon[curr_time][num]
                i, j, k = int(xBig/xStepGrid), int(yBig/yStepGrid), int(zBig/zStepGrid)
                speeds = getSpeedProjection(positionCarbon[curr_time][num])

                for l in range(positionCarbon.shape[2]):
                    positionCarbon[p_next_time(time)][num][l] = positionCarbon[curr_time][num][l]
                for dim in range(3):
                    v = speeds[dim]
                    r = positionCarbon[curr_time][num][dim]/Ml
                    E = tensionCorpusculActingCarbon[num][dim]/Mee
                    delta, Z, epsilonC = 1.0, 1.0, 1.0

                    y0 = [r, v, E, delta, Z, epsilonC]
                    res = odeint(f, y0, t)
                    positionCarbon[p_next_time(time)][num][dim] = res[-1][0]*Ml
                    positionCarbon[p_next_time(time)][num][5+dim] = res[-1][1]
            for num in range(positionElectron.shape[1]):
                for l in range(positionElectron.shape[2]):
                    positionElectron[p_next_time(time)][num][l] = positionElectron[curr_time][num][l]
            for num in range(positionHelium.shape[1]):
                for l in range(positionHelium.shape[2]):
                    positionHelium[p_next_time(time)][num][l] = positionHelium[curr_time][num][l]
            print('THIS IS END OF RKF')

            time += 1
            print('time = {} '.format(time))

            directory = make_dir(prefix, 'by_time')    
            n = (xDimensionGrid, yDimensionGrid, zDimensionGrid)
            fig = plt.figure()
            ax = fig.add_subplot(211, projection='3d')
            ax2 = fig.add_subplot(212)
            plt.title('time = {} '.format(time))
            ax = make_3d_plot_with_speed(ax, positionCarbon, p_prev_time(time), (xStepGrid, yStepGrid, zStepGrid), n, plot_type='COORD')
            ax2 = make_2d_plot_with_tracks(ax2, positionCarbon, p_prev_time(time), (xStepGrid, yStepGrid, zStepGrid), n)
            # plt.show()
            # plt.savefig("./picts/carbon_3d_speed_{0:04d}.png".format(time))
            plt.savefig("{}/carbon_two_coord_{:04d}.png".format(directory, time))
            plt.clf()
            # fig = None

            # directory = make_dir(prefix, 'inten')    
            # fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
            # for i, ax in zip(inten, [ax1, ax2, ax3]):
            #     ax.plot([k for k,_,_ in i], color='red')
            # plt.savefig("{}/int_time={:04d}".format(directory, time))
            # plt.clf()


            # directory = make_dir(prefix, 'phi')    
            # n = next_phi.shape;
            # plot_phi = np.zeros((n[0], n[1]))
            # for k in range(n[2]):
            #     for i in range(n[0]):
            #         for j in range(n[1]):
            #             plot_phi[i][j] = next_phi[i][j][k]

            #     plt.contourf(plot_phi.T, cmap=plt.cm.flag)
            #     plt.colorbar()
            #     plt.savefig("{}/phi_time={:04d}_z={:02d}".format(directory, time, k))
            #     plt.clf()

            # directory = make_dir(prefix, 'intensity')    
            # n = intensity.shape;
            # plot_intensity = np.zeros((n[0], n[1]))
            # for k in range(n[2]):
            #     for i in range(n[0]):
            #         for j in range(n[1]):
            #             # plot_intensity[i][j] = np.sqrt(np.sum([intensity[i][j][k][l]**2 for l in range(3)]))
            #             plot_intensity[i][j] = intensity[i][j][k][0]

            #     plt.contourf(plot_intensity.T, cmap=plt.cm.flag)
            #     plt.colorbar()
            #     plt.savefig("{}/intensity_time={:04d}_z={:02d}".format(directory, time, k))
            #     plt.clf()

            # directory = make_dir(prefix, 'tension')    
            # tension = tensionCorpusculActingCarbon
            # n = tension.shape;
            # plot_tension = np.zeros(n[0])
            # plot_tension_x = np.zeros(n[0])
            # plot_tension_y = np.zeros(n[0])
            # plot_tension_z = np.zeros(n[0])
            # for i in range(n[0]):
            #     plot_tension[i] = np.sqrt(np.sum([tension[i][l]**2 for l in range(3)]))
            #     plot_tension_x[i] = tension[i][0]
            #     plot_tension_y[i] = tension[i][1]
            #     plot_tension_z[i] = tension[i][2]

            # fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
            # ax1.plot(plot_tension_x)
            # ax2.plot(plot_tension_y)
            # ax3.plot(plot_tension_z)
            # plt.savefig("{}/tension_by_x_y_z_time={:04d}".format(directory, time))
            # plt.clf()

            # plt.plot(plot_tension)
            # plt.savefig("{}/tension_by_sum(x,y,z)_time={:04d}".format(directory, time))
            # plt.clf()


            
            curr_time = p_prev_time(time)
            for p, l in zip(listen_particles, range(len(listen_particles))):
                plot_data[l] += [positionCarbon[curr_time][p].copy()]
            # print(plot_data)
    except IndexError:
        print('IndexError')
        pass

    directory = make_dir(prefix, 'graphs')
    colors = ['red', 'green', 'blue']
    for pd, p in zip(plot_data, listen_particles):    
        fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
        ax1.plot( [i[0] for i in pd])#, color=colors[p])
        ax1.set_title("x from time particle num={}".format(p))
        ax2.plot( [i[1] for i in pd])#, color=colors[p])
        ax2.set_title("y from time particle num={}".format(p))
        ax3.plot( [i[2] for i in pd])#, color=colors[p])
        ax3.set_title("z from time particle num={}".format(p))
        plt.savefig("{}/positions_50_50_pn={}".format(directory, p))
        plt.clf()

        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)
        ax1.plot([i[0] for i in pd], [i[5]*Mnue if np.abs(i[5]*Mnue)>1E-30 else 0 for i in pd])#, color=colors[p])
        ax1.set_title("v_x from x particle num={}".format(p))
        ax2.plot([i[0] for i in pd], [i[6]*Mnue if np.abs(i[6]*Mnue)>1E-30 else 0 for i in pd])#, color=colors[p])
        ax2.set_title("v_y from x particle num={}".format(p))
        ax3.plot([i[0] for i in pd], [i[7]*Mnue if np.abs(i[7]*Mnue)>1E-30 else 0 for i in pd])#, color=colors[p])
        ax3.set_title("v_z from x particle num={}".format(p))
        ax4.plot([i[0] for i in pd], [np.sqrt(i[5]**2 + i[6]**2 + i[7]**2)*Mnue for i in pd])#, color=colors[p])
        ax4.set_title("|v| from x particle num={}".format(p))
        plt.savefig("{}/speeds_50_50_pn={}.png".format(directory, p))
        plt.clf()

        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)
        ax1.plot([i for i in range(len(pd))], [i[5]*Mnue if np.abs(i[5]*Mnue)>1E-30 else 0 for i in pd])#, color=colors[p])
        ax1.set_title("v_x from time particle num={}".format(p))
        ax2.plot([i for i in range(len(pd))], [i[6]*Mnue if np.abs(i[6]*Mnue)>1E-30 else 0 for i in pd])#, color=colors[p])
        ax2.set_title("v_y from time particle num={}".format(p))
        ax3.plot([i for i in range(len(pd))], [i[7]*Mnue if np.abs(i[7]*Mnue)>1E-30 else 0 for i in pd])#, color=colors[p])
        ax3.set_title("v_z from time particle num={}".format(p))
        ax4.plot([i for i in range(len(pd))], [np.sqrt(i[5]**2 + i[6]**2 + i[7]**2)*Mnue for i in pd])#, color=colors[p])
        ax4.set_title("|v| from time particle num={}".format(p))
        plt.savefig("{}/speeds_50_50_from_time_pn={}.png".format(directory, p))
        plt.clf()

        fig, (ax1, ax2) = plt.subplots(2, sharex=True)
        ax1.plot([i[0] for i in pd], [i[1] for i in pd])
        ax1.set_title("y from x particle num={}".format(p))
        ax2.plot([i[0] for i in pd], [i[2] for i in pd])#, color=colors[p])
        ax2.set_title("z from x particle num={}".format(p))
        plt.savefig("{}/positions(position)_50_50_pn={}".format(directory, p))

        fig, (ax1) = plt.subplots(1)
        ax1.plot([i[1] for i in pd], [i[2] for i in pd])#, color=colors[p])
        ax1.set_title("z from y particle num={}".format(p))
        plt.savefig("{}/positions(position)_2_50_50_pn={}".format(directory, p))
        plt.clf()

import os
import time
import shutil

def make_dir_prefix():
    prefix = './picts/picts_cuda_{}'.format(time.strftime('%Y%m%d%H%M%S'))
    # prefix = './picts_{}'.format(time.strftime('%Y%m%d%H'))#%M%S'))
    try:
        os.mkdir(prefix)
    except:
        pass
    shutil.copy('./constant.py', prefix)
    return prefix

def make_dir(prefix, dir_name):
    path = os.path.join(prefix, dir_name)
    try:
        os.mkdir(path)
    except:
        pass
    return path




if __name__ == '__main__':
    prefix = make_dir_prefix()
    # print(make_dir(prefix, 'fff'))
    main(prefix)
