# -*- coding: utf-8 -*-

import os
import time
import shutil
import scipy
import scipy.stats
from scipy.integrate import odeint
import numpy as np
from numpy.random import rand

from constant import *
from plot import *
from establishing_method import *
from establishing_method import potential_establish_method_cuda as em
from collision import *
from dump import *

__author__ = 'fedul0x'  

def p_time(absolute_time):
    ''' pseuso time '''
    return absolute_time % (DATA_IN_MEMORY_TIME)


def p_next_time(absolute_time):
    return (absolute_time + 1) % (DATA_IN_MEMORY_TIME)


def p_prev_time(absolute_time):
    return (absolute_time - 1) % (DATA_IN_MEMORY_TIME)


def f(y, t):
    delta, Z, epsilonC = 1.0, 1.0, 1.0
    r, v, E = y[0], y[1], y[2]
    f = [np.sqrt(delta)*v, np.sqrt(delta)*Z/2/epsilonC*E, 0]
    return f


class MWranomizer(scipy.stats.rv_continuous):

    def __init__(self, n, m, T):
        ''' Maxwell idstribution randomizer '''
        scipy.stats.rv_continuous.__init__(self)
        self.a = 0
        self.n, self.m, self.T = n, m, T

    def _pdf(self, v):
        k = BOLTZMANN_CONSTANT
        n, m, T = self.n, self.m, self.T
        return n*(m/(2*np.pi*k*T))**(3/2)*np.exp(-m*v*v/(2*k*T))*4*np.pi*v**2

class MWranomizer2:

    def __init__(self, n, m, T):
        ''' Maxwell idstribution randomizer '''
        self.a = 0
        self.n, self.m, self.T = n, m, T

    def rvs(self, size=1):
        return scipy.stats.maxwell.rvs(scale=CARBONS_MAX_SPEED/1.5957691216057308, size=size)

def get_component(particle, b=0, n=3):
    return [particle[i] for i in range(b, b+n)]

def get_collision_energy(part1, part2):
    xc, yc, zc, radiusc, chrg, vcx, vcy, vcz, _ = part1
    xc2, yc2, zc2, radiusc2, chrg2, vcx2, vcy2, vcz2, _ = part2
    m1 = CARBONS_MASS*radiusc**3 / PARTICLE_COHESION / CARBONS_RADIUS**3
    m2 = CARBONS_MASS*radiusc2**3 / PARTICLE_COHESION / CARBONS_RADIUS**3
    v1 = CARBON_SPEED_DIMENSIONLESS_UNIT*np.sqrt(vcx**2 + vcy**2 + vcz**2)
    v2 = CARBON_SPEED_DIMENSIONLESS_UNIT*np.sqrt(vcx2**2 + vcy2**2 + vcz2**2)
    e1 = m1 * v1**2 / 2 / CARBONS_NUMBER * AVOGADRO_CONSTANT
    e2 = m2 * v2**2 / 2 / CARBONS_NUMBER * AVOGADRO_CONSTANT
    return (e1, e2)
# def getSpeedProjection(particle):
#     # speed = particle[5]
#     # devider = np.sqrt(3)
#     # return (speed*Kvu￼1/devider, speed*Kvu2/devider, speed*Kvu2/devider)
#     return (particle[5], particle[6], particle[7])


def spread_position(center, number):
    # x, y, z, = center
    # x_big, y_big, z_big = 0, 0, 0
    # for _ in range(int(number)):
    #     x_big += x + rand() * X_STEP
    #     y_big += y + rand() * Y_STEP
    #     z_big += z + rand() * Z_STEP
    # return x_big/number, y_big/number, z_big/number
    x, y, z, = center
    big_pos = []
    for p, step in zip(center, (X_STEP, Y_STEP, Z_STEP,)):
        big_pos += [np.random.normal(loc=p + step/2, scale=(step/6)**2)]
    return big_pos

def spread_carbon_position(center, number):
    # x, y, z, = center
    # x_big, y_big, z_big = 0, 0, 0
    # for _ in range(int(number)):
    #     x_big += x + rand() * X_STEP
    #     y_big += y + rand() * Y_STEP
    #     z_big += z + rand() * Z_STEP
    # return x_big/number, y_big/number, z_big/number
    x, y, z, = center
    p, step = center[0], X_STEP*X_STEP_NUMBER_GRID/4
    big_pos = [np.abs(np.random.normal(loc=p + step/2, scale=(step*30)**2))]
    # print('carbon x = {} < {}*25 = {}'.format(big_pos[0], X_STEP, X_STEP**X_STEP_NUMBER_GRID/4))
    for p, step in zip(center[1:], (Y_STEP, Z_STEP,)):
        big_pos += [np.random.normal(loc=p + step/2, scale=(step/6)**2)]
    return big_pos

def spread_speed(randomizer, dimensionless=1):
    speed = randomizer.rvs(size=1)[0]/dimensionless
    mu, sigma = 0, SPEED_DISTRIBUTION_ANGLE_SIGMA
    angle = np.random.normal(mu, sigma, 2)
    angle[1] = angle[1] + np.pi/2
    x_speed = speed*np.cos(angle[0])
    y_speed = speed*np.cos(angle[1])
    z_speed = np.abs((x_speed/np.cos(angle[0]))**2 - x_speed**2 - y_speed**2)
    return x_speed, y_speed, np.sqrt(z_speed)


def spread_charge(center, big_center, charge):
    x, y, z, = center
    i, j, k = 0, 0, 0
    x_big, y_big, z_big = big_center
    v = X_STEP * Y_STEP * Z_STEP
    x_neighbor, y_neighbor, z_neighbor = x_big - x, y_big - y, z_big - z
    x_opposite, y_opposite, z_opposite = \
        X_STEP - x_big + x, Y_STEP - y_big + y, Z_STEP - z_big + z

    res = []
    res += [(i, j, k, charge * x_opposite * y_opposite * z_opposite / v)]
    res += [(i + 1, j, k, charge * x_neighbor * y_opposite * z_opposite / v)]
    res += [(i, j + 1, k, charge * x_opposite * y_neighbor * z_opposite / v)]
    res += [(i, j, k + 1, charge * x_opposite * y_opposite * z_neighbor / v)]
    res += [(i + 1, j + 1, k, charge * x_neighbor * y_neighbor * z_opposite / v)]
    res += [(i, j + 1, k + 1, charge * x_opposite * y_neighbor * z_neighbor / v)]
    res += [(i + 1, j, k + 1, charge * x_neighbor * y_opposite * z_neighbor / v)]
    res += [(i + 1, j + 1, k + 1, charge * x_neighbor * y_neighbor * z_neighbor / v)]

    return res


def spread_tension(center, big_center, intensity):
    x, y, z, = center
    x_big, y_big, z_big = big_center
    m = intensity
    i, j, k = int(x_big/X_STEP), int(y_big/Y_STEP), int(z_big/Z_STEP)
    v = X_STEP * Y_STEP * Z_STEP
    X_neighbor, Y_neighbor, Z_neighbor = x_big - x, y_big - y, z_big - z
    x_opposite, y_opposite, z_opposite = \
        X_STEP - x_big + x, Y_STEP - y_big + y, Z_STEP - z_big + z
    tension = [0, 0, 0]
    for index in range(3):
        tension[index] = \
            m[i][j][k][index] * x_opposite * y_opposite * z_opposite / v +\
            m[i+1][j][k][index] * X_neighbor * y_opposite * z_opposite / v +\
            m[i][j+1][k][index] * x_opposite * Y_neighbor * z_opposite / v +\
            m[i][j][k+1][index] * x_opposite * y_opposite * Z_neighbor / v +\
            m[i+1][j+1][k][index] * X_neighbor * Y_neighbor * z_opposite / v +\
            m[i][j+1][k+1][index] * x_opposite * Y_neighbor * Z_neighbor / v +\
            m[i+1][j][k+1][index] * X_neighbor * y_opposite * Z_neighbor / v +\
            m[i+1][j+1][k+1][index] * X_neighbor * Y_neighbor * Z_neighbor / v

    return tension


def dif(x, y, z, step):
    return (x - 2 * y + z) / step ** 2


def main():
    # Границы сетки
    x_range = np.linspace(0, X_DIMENSION_GRID, X_STEP_NUMBER_GRID + 1)
    y_range = np.linspace(0, Y_DIMENSION_GRID, Y_STEP_NUMBER_GRID + 1)
    z_range = np.linspace(0, Z_DIMENSION_GRID, Z_STEP_NUMBER_GRID + 1)
    n = (x_range.shape[0]-1) * (y_range.shape[0]-1) * (z_range.shape[0]-1)
    cpl = int(n/(x_range.shape[0] - 1))
    # time, i, x, y, z, radius, charge, speedx, speedy, speedz, guid
    electron = np.empty([DATA_IN_MEMORY_TIME, n, 6+2+1])
    carbon = \
        np.empty([DATA_IN_MEMORY_TIME, CARBON_LAYERS_NUMBER*cpl, 6+2+1])
    helium = np.empty([DATA_IN_MEMORY_TIME, n, 6+2+1])
    # Сразу крупные частицы
    num = 0  # particle number
    carbon_num = 0  # carbon number
    carbons_in_process = cpl  # carbon particles across layers
    absolute_time = 0  # time scale position
    # TODO оптимизировать получения случвела Максвелла

    start = time.time()
    csdu, esdu, hsdu = CARBON_SPEED_DIMENSIONLESS_UNIT, \
        ELECTRON_SPEED_DIMENSIONLESS_UNIT, HELLIUM_SPEED_DIMENSIONLESS_UNIT
    big_carbon_radius = \
        (PARTICLE_COHESION * CARBONS_RADIUS ** 3 * CARBONS_NUMBER) ** (1.0 / 3.0)
    big_carbon_charge = CARBONS_CHARGE * CARBONS_NUMBER
    big_electron_radius = \
        (PARTICLE_COHESION * ELECTRONS_RADIUS ** 3 * ELECTRONS_NUMBER)**(1.0 / 3.0)
    big_electron_charge = ELECTRONS_CHARGE * ELECTRONS_NUMBER
    big_helium_radius = \
        (PARTICLE_COHESION * HELIUMS_RADIUS ** 3 * HELIUMS_NUMBER) ** (1.0 / 3.0)
    big_helium_charge = HELIUMS_CHARGE * HELIUMS_NUMBER
    prev_phi, next_phi = [], []
    # Restore from dump particles and potential field
    is_dumped = True
    data = restore_from_dump(DUMP_FOLDER, CONSTANT_VALUES)
    if not (data is None):
        print('Resore from dump')
        carbon, electron, helium, prev_phi, next_phi = data
        print('Particle distribution extraction from dump elapsed time = {}'.format(time.time()-start))
    else:
        print('Particle distribution generation')
        electron_randomizer = MWranomizer(n=1, m=ELECTRONS_MASS, T=TEMPERATURE)
        carbon_randomizer = MWranomizer(n=1, m=CARBONS_MASS, T=TEMPERATURE)
        helium_randomizer = MWranomizer(n=1, m=HELIUMS_MASS, T=HELIUMS_TEMPERATURE)
        # Carbon distribution
        is_dumped = False
        for _ in range(CARBON_LAYERS_NUMBER):
            for y, j in zip(y_range[:-1], range(y_range.shape[0] - 1)):
                for z, k in zip(z_range[:-1], range(z_range.shape[0] - 1)):
                    cell = (x_range[:-1][0]*25, y, z)
                    x_big, y_big, z_big = \
                        spread_carbon_position(cell, CARBONS_NUMBER)
                    x_speed, y_speed, z_speed = \
                        spread_speed(carbon_randomizer, dimensionless=csdu)
                    for v, l in zip([x_big, y_big, z_big, big_carbon_radius, big_carbon_charge, x_speed, y_speed, z_speed, carbon_num], range(carbon.shape[2])):
                        carbon[absolute_time][carbon_num][l] = v
                    carbon_num += 1
        # Electron and helium distribution
        for y, j in zip(y_range[:-1], range(y_range.shape[0] - 1)):
            for z, k in zip(z_range[:-1], range(z_range.shape[0] - 1)):
                for x, i in zip(x_range[:-1], range(x_range.shape[0] - 1)):
                    cell = (x, y, z)
                    x_big, y_big, z_big = \
                        spread_position(cell, ELECTRONS_NUMBER)
                    x_speed, y_speed, z_speed = \
                        spread_speed(electron_randomizer, dimensionless=esdu)
                    for v, l in zip([x_big, y_big, z_big, big_electron_radius, big_electron_charge, x_speed, y_speed, z_speed, num], range(electron.shape[2])):
                        electron[absolute_time][num][l] = v
                    x_big, y_big, z_big = \
                        spread_position(cell, HELIUMS_NUMBER)
                    x_speed, y_speed, z_speed = \
                        spread_speed(helium_randomizer, dimensionless=hsdu)
                    for v, l in zip([x_big, y_big, z_big, big_helium_radius, big_helium_charge, x_speed, y_speed, z_speed, num], range(helium.shape[2])):
                        helium[absolute_time][num][l] = v
                    num += 1
        print('Particle distribution elapsed time = {}'.format(time.time()-start))

    # MODELING CYCLE BEGIN
    num = 0  # number of particle
    absolute_time = 0  # time scale position
    max_carbon_guid = carbon.shape[1]
    marshal = lambda x: '{} = {}'.format(x[0], x[1])
    constkeys = sorted(CONSTANT_VALUES.items())
    constants = map(marshal, constkeys)
    db_log = DbConnection(DB_FILE, '\n'.join(constants))
    try:
        while (absolute_time < MODELING_TIME) and (carbons_in_process > 0):
            curr_time = p_time(absolute_time)
            # Charge in nodes calculations
            start = time.time()
            sizes = [x_range.shape[0], y_range.shape[0], z_range.shape[0]]
            electron_charge_grid = np.zeros(sizes)
            carbon_charge_grid = np.zeros(sizes)
            helium_charge_grid = np.zeros(sizes)
            positions = [electron, carbon, helium]
            sizes = [electron.shape[1], carbons_in_process, helium.shape[1]]
            grids = [electron_charge_grid, carbon_charge_grid, helium_charge_grid]
            names = ['electron', 'carbon', 'helium']
            carbon_final, tension_final = [], []
            for grid, position, name, size in zip(grids, positions, names, sizes):
                n = grid.shape
                for num in range(size):
                    x_big, y_big, z_big, _, charge, _, _, _, _ = position[curr_time][num]
                    try:
                        i, j, k = \
                            int(x_big/X_STEP), int(y_big/Y_STEP), int(z_big/Z_STEP)
                    except Exception:
                        print('{}/{}, {}/{}, {}/{}'.format(x_big, X_STEP, y_big, Y_STEP, z_big, Z_STEP))
                    i, j, k = \
                        int(x_big/X_STEP), int(y_big/Y_STEP), int(z_big/Z_STEP)
                    # Carbon logging when it reaches the end of the simulation area
                    if name == 'carbon':
                        if (i<0) or (j<0) or (k<0) or (i>n[0]-2) or (j>n[1]-2) or (k>n[2]-2):
                            continue
                    x, y, z = \
                        i*X_STEP, j*Y_STEP, k*Z_STEP
                    patch = spread_charge((x, y, z), (x_big, y_big, z_big), charge)
                    for p in patch:
                        grid[i+p[0]][j+p[1]][k+p[2]] += p[3]
            end = time.time()
            print('Calcing charge in cells elapsed time = {}'.format(end-start))

            # The boundary conditions for the calculation of the potential
            start = time.time()
            n = (x_range.shape[0]+2, y_range.shape[0]+2, z_range.shape[0]+2)
            phi = POTENTIAL_BOUND_VALUE
            ecg, ccg, hcg = \
                electron_charge_grid, carbon_charge_grid, helium_charge_grid
            # TODO easying
            start2 = time.time()
            if FAST_ESTABLISHING_METHOD:
                if prev_phi == [] and next_phi == []:
                    prev_phi, next_phi, ro = \
                        make_boundary_conditions(phi, n, ecg, ccg, hcg)
                else:
                    _, _, ro = \
                        make_boundary_conditions(phi, n, ecg, ccg, hcg)
            else:
                prev_phi, next_phi, ro = \
                    make_boundary_conditions(phi, n, ecg, ccg, hcg)
            end2 = time.time()
            print('Preparation elapsed time = {}'.format(end2-start2))
            
            # Potential establish method
            ema = ESTABLISHING_METHOD_ACCURACY
            prev_phi, next_phi = \
                em(prev_phi, next_phi, ro, epsilon=ema)
            # make_potential_plot_file('./picts/', next_phi, curr_time)
            # Dump particles and potential field
            if not is_dumped:
                data = (carbon, electron, helium, prev_phi, next_phi)
                save_to_dump(DUMP_FOLDER, CONSTANT_VALUES, data)
                is_dumped = True
            end = time.time()
            # make_potential_plot_file('./picts/', next_phi, absolute_time)
            print('Establishing method elapsed time = {}'.format(end-start))
            start = time.time()
            n = next_phi.shape
            intensity = np.zeros([n[0]-2, n[1]-2, n[2]-2, 3])
            inten = [[], [], []]
            for i in range(1, n[0]-1):
                for j in range(1, n[1]-1):
                    for k in range(1, n[2]-1):
                        intensity[i-1][j-1][k-1][0] = (next_phi[i - 1][j][k] - next_phi[i + 1][j][k]) / 2 / X_STEP
                        intensity[i-1][j-1][k-1][1] = (next_phi[i][j - 1][k] - next_phi[i][j + 1][k]) / 2 / Y_STEP
                        intensity[i-1][j-1][k-1][2] = (next_phi[i][j][k - 1] - next_phi[i][j][k + 1]) / 2 / Z_STEP
                        inten[0] += [(intensity[i-1][j-1][k-1][0], next_phi[i - 1][j][k], next_phi[i + 1][j][k])]
                        inten[1] += [(intensity[i-1][j-1][k-1][1], next_phi[i][j - 1][k], next_phi[i][j + 1][k])]
                        inten[2] += [(intensity[i-1][j-1][k-1] [2], next_phi[i][j][k - 1], next_phi[i][j][k + 1])]
            end = time.time()
            # make_inten_plot_file('./picts/', inten, absolute_time)
            print('Intensity calcing elapsed time = {}'.format(end-start))
            # Calcing tension by each particle
            start = time.time()
            electron_tension = np.empty([electron.shape[1], 3])
            carbon_tension = np.empty([carbons_in_process, 3])
            helium_tension = np.empty([helium.shape[1], 3])
            particles = [electron, carbon, helium]
            sizes = [electron.shape[1], carbons_in_process, helium.shape[1]]
            tensions = [electron_tension, carbon_tension, helium_tension]
            for position, p, size in zip(particles, range(len(tensions)), sizes):
                for num in range(size):
                    x_big, y_big, z_big = \
                        get_component(position[curr_time][num])
                    i, j, k = \
                        int(x_big/X_STEP), int(y_big/Y_STEP), int(z_big/Z_STEP)
                    x, y, z = i*X_STEP, j*Y_STEP, k*Z_STEP

                    # print('x_big, y_big, z_big = {}, {}, {}'.format(x_big, y_big, z_big))
                    # print('i, j, k = {}, {}, {}'.format(i, j, k))
                    # print('x, y, z = {}, {}, {}'.format(x, y, z))
                    try:
                        tension = spread_tension((x, y, z), (x_big, y_big, z_big), intensity)
                    except:
                        print(num)
                        return
                    for t in range(3):
                        tensions[p][num][t] = tension[t]
            # make_tension_plot_file('./picts/', carbon_tension, absolute_time)
            end = time.time()
            print('Tension by particle calcing elapsed time = {}'.format(end-start))
            
            # Solutions of differential equations
            start = time.time()
            t = np.linspace(0, TIME_STEP / Mt, 2)
            curr_time = p_time(absolute_time)
            for num in range(carbons_in_process):
                x_big, y_big, z_big = \
                    get_component(carbon[curr_time][num])
                i, j, k = \
                    int(x_big/X_STEP), int(y_big/Y_STEP), int(z_big/Z_STEP)
                speeds = get_component(carbon[curr_time][num], b=5)
                for l in range(carbon.shape[2]):
                    carbon[p_next_time(absolute_time)][num][l] = carbon[curr_time][num][l]
                for dim in range(3):
                    # import pdb
                    # pdb.set_trace()
                    v = speeds[dim]
                    r = carbon[curr_time][num][dim]/SPACE_DIMENSIONLESS_UNIT
                    E = carbon_tension[num][dim]/INTENSITY_DIMENSIONLESS_UNIT
                    y0 = [r, v, E]
                    fff = True
                    while fff:
                        res = odeint(f, y0, t)
                        fff = (res[-1][0] == np.NaN) or (res[-1][1] == np.NaN)
                    # if num == 0:
                    #     print(res)
                    if res[-1][0] != np.NaN:
                        carbon[p_next_time(absolute_time)][num][dim] = res[-1][0]*SPACE_DIMENSIONLESS_UNIT
                    else:
                        carbon[p_next_time(absolute_time)][num][dim] = carbon[curr_time(absolute_time)][num][dim]
                    if res[-1][1] != np.NaN:
                        carbon[p_next_time(absolute_time)][num][5+dim] = res[-1][1]
                    else:
                        carbon[p_next_time(absolute_time)][num][5+dim] = carbon[curr_time(absolute_time)][num][5+dim]
            end = time.time()
            print('Differential equations solution elapsed time = {}'.format(end-start))
            
# Copy another particles arrays
            for num in range(carbons_in_process, carbon.shape[1]):
                for l in range(carbon.shape[2]):
                    carbon[p_next_time(absolute_time)][num][l] = \
                        carbon[curr_time][num][l]
            for num in range(electron.shape[1]):
                for l in range(electron.shape[2]):
                    electron[p_next_time(absolute_time)][num][l] = \
                        electron[curr_time][num][l]
            for num in range(helium.shape[1]):
                for l in range(helium.shape[2]):
                    helium[p_next_time(absolute_time)][num][l] = \
                        helium[curr_time][num][l]
                        
            # positions = [electron, carbon, helium]
            offsets, carbon_final, tension_final = [], [], []
            curr_time = p_next_time(absolute_time)
            n = carbon_charge_grid.shape
            # print(n)
            # return
            for num in range(carbons_in_process):
                x_big, y_big, z_big, _, charge, _, _, _, _ = carbon[curr_time][num]
                i, j, k = \
                    int(x_big/X_STEP), int(y_big/Y_STEP), int(z_big/Z_STEP)
                # Carbon logging when it reaches the end of the simulation area
                # if num == 23:
                #     print('23: i, j, k = {}, {}, {}'.format(i, j, k))
                if (i<0) or (j<0) or (k<0) or (i>n[0]-2) or (j>n[1]-2) or (k>n[2]-2):
                    offsets += [num]
                    carbon_final += [carbon[curr_time][num]]
                    tension_final += [carbon_tension[num]]
            
            end = time.time()
            print('Calcing charge in cells elapsed time = {}'.format(end-start))




            # Collision finding
            start = time.time()
            сarbon_collisions, electron_collisions, helium_collisions = [], [], []
            collision_offsets = []
            if absolute_time != 0:
                prev_time = p_time(absolute_time)
                curr_time = p_next_time(absolute_time)
                сarbon_collisions = []
                сarbon_collisions_tmp = find_carbon_collision_rust(carbon, carbons_in_process, curr_time, prev_time)
                for col in сarbon_collisions_tmp:
                    if col[0] in offsets or col[1] in offsets:
                        continue
                    e1, e2 = get_collision_energy(carbon[curr_time][col[0]], carbon[curr_time][col[1]])
                    xc, yc, zc, radiusc, chrg, vcx, vcy, vcz, _ = carbon[curr_time][col[0]]
                    xc2, yc2, zc2, radiusc2, chrg2, vcx2, vcy2, vcz2, _ = carbon[curr_time][col[1]]
                    m1 = CARBONS_MASS*radiusc**3 / PARTICLE_COHESION / CARBONS_RADIUS**3
                    m2 = CARBONS_MASS*radiusc2**3 / PARTICLE_COHESION / CARBONS_RADIUS**3
                    if e1 + e2 > 348000:
                        import matplotlib.pyplot as plt
                        pxc, pyc, _, _, _, _, _, _, _ = carbon[prev_time][col[0]]
                        pxc2, pyc2, _, _, _, _, _, _, _ = carbon[prev_time][col[1]]
                        # plt.plot([pxc, xc], [pyc, yc], '-b')
                        # plt.plot([pxc2, xc2], [pyc2, yc2], '-r')
                        r = (((m1 + m2) * PARTICLE_COHESION * CARBONS_RADIUS**3) / CARBONS_MASS)**(1.0/3)
                        carbon[curr_time][col[0]][3] = r
                        carbon[curr_time][col[0]][4] = chrg + chrg2
                        carbon[curr_time][col[0]][5] = vcx + vcx2
                        carbon[curr_time][col[0]][6] = vcy + vcy2
                        carbon[curr_time][col[0]][7] = vcz + vcz2
                        carbon[curr_time][col[0]][-1] = max_carbon_guid
                        max_carbon_guid += 1
                        # plt.plot([carbon[curr_time][col[0]][0]], [carbon[curr_time][col[0]][1]], 'or')
                        # plt.xlim(0, 0.001)
                        # plt.ylim(0, 0.01)
                        # plt.savefig('./picts_4/{}_{}.png'.format(col[0], col[1]))
                        # plt.close()
                        # fff = open('./collision.txt', 'a')
                        # fff.write('{}\n'.format(xc))
                        # fff.close()
                        collision_offsets += [col[1]]
                        # continue
                    else:
                        carbon[curr_time][col[0]][5] = vcx2
                        carbon[curr_time][col[0]][6] = vcy2
                        carbon[curr_time][col[0]][7] = vcz2
                        carbon[curr_time][col[1]][5] = vcx
                        carbon[curr_time][col[1]][6] = vcy
                        carbon[curr_time][col[1]][7] = vcz
                    сarbon_collisions += [col]

            end = time.time()
            curr_time = p_time(absolute_time)
            print('Collisions finding elapsed time = {}'.format(end-start))

            curr_time = p_next_time(absolute_time)
            def enrg(x):
                p1 = carbon[curr_time][x[0]]
                p2 = carbon[curr_time][x[1]]
                e = get_collision_energy(p1, p2)
                return (x[0], x[1], sum(e))
            сarbon_collisions = map(enrg, сarbon_collisions)
            db_log.new_iteration(absolute_time*TIME_STEP)
            db_log.new_final(carbon_final, tension_final)
            # db_log.new_collision(сarbon_collisions)
            db_log.new_collision_2(carbon, сarbon_collisions, p_next_time(absolute_time))
            curr_time = p_time(absolute_time)


            curr_time = p_next_time(absolute_time)
            if offsets:
                v1, v2 = min(offsets), min(offsets)
                offset, i = 0, v1
                while i < carbon.shape[1]:
                    while i+offset in offsets:
                        offset += 1
                    if i + offset < carbon.shape[1]:
                        carbon[curr_time][i] = carbon[curr_time][i+offset]
                    i += 1
                carbons_in_process -= len(offsets)
                print('offsets: {}'.format(len(offsets)))
            curr_time = p_time(absolute_time)

            curr_time = p_next_time(absolute_time)
            if collision_offsets:
                v1, v2 = min(collision_offsets), min(collision_offsets)
                offset, i = 0, v1
                while i < carbon.shape[1]:
                    while i+offset in collision_offsets:
                        offset += 1
                    if i + offset < carbon.shape[1]:
                        carbon[curr_time][i] = carbon[curr_time][i+offset]
                        # carbon[curr_time][i][-1] = max_carbon_guid
                        # max_carbon_guid += 1
                    i += 1
                carbons_in_process -= len(collision_offsets)
                print('collision_offsets: {}'.format(len(collision_offsets)))

            db_log.new_particle_2(carbon, p_next_time(absolute_time), carbons_in_process)
            curr_time = p_time(absolute_time)

            




            absolute_time += 1
            print('PREV carbons_in_process = {}'.format(carbons_in_process))
            if absolute_time < CARBON_LAYERS_NUMBER:
                carbons_in_process += cpl
            print('absolute time = {} '.format(absolute_time))
            print('carbons_in_process = {}'.format(carbons_in_process))
            # print(carbon[p_prev_time(absolute_time)][23])
            # print(carbon[p_time(absolute_time)][23])
            # print(offsets)
            # print(collision_offsets)
    # except IndexError:
    #     print('IndexError')
    #     pass
    except KeyboardInterrupt:
        print('KeyboardInterrupt')
        pass
    # except OSError:
    #     print('OSError')
    #     pass


if __name__ == '__main__':
    main()
            # Calculation of tension