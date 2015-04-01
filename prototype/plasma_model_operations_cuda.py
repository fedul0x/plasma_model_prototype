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


def p_time(time):
    ''' pseuso time '''
    return time % (DATA_IN_MEMORY_TIME)


def p_next_time(time):
    return (time + 1) % (DATA_IN_MEMORY_TIME)


def p_prev_time(time):
    return (time - 1) % (DATA_IN_MEMORY_TIME)


def f(y, t):
    delta, Z, epsilonC = 1.0, 1.0, 1.0
    r, v, E = y[0], y[1], y[2]
    f = [np.sqrt(delta)*v, np.sqrt(delta)*Z/2/epsilonC*E]
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


def get_component(particle, b=0, n=3):
    return [particle[i] for i in range(b, b+n)]


# def getSpeedProjection(particle):
#     # speed = particle[5]
#     # devider = np.sqrt(3)
#     # return (speed*Kvu￼1/devider, speed*Kvu2/devider, speed*Kvu2/devider)
#     return (particle[5], particle[6], particle[7])


def spread_position(center, number):
    x, y, z, = center
    x_big, y_big, z_big = 0, 0, 0
    for _ in range(number):
        x_big += x + rand() * X_STEP
        y_big += y + rand() * Y_STEP
        z_big += z + rand() * Z_STEP
    return x_big/number, y_big/number, z_big/number

def spread_speed(randomizer, dimensionless=1):
    speed = randomizer.rvs(size=1)[0]/dimensionless
    # x_speed, y_speed, z_speed = speed, 0, 0
    mu, sigma = 0, np.pi/30
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


def main(prefix):
    # Границы сетки
    x_range = np.linspace(0, X_DIMENSION_GRID, X_STEP_NUMBER_GRID + 1)
    y_range = np.linspace(0, Y_DIMENSION_GRID, Y_STEP_NUMBER_GRID + 1)
    z_range = np.linspace(0, Z_DIMENSION_GRID, Z_STEP_NUMBER_GRID + 1)
    n = (x_range.shape[0]-1) * (y_range.shape[0]-1) * (z_range.shape[0]-1)
    # time, i, x, y, z, radius, charge, speedx, speedy, speedz
    positionElectron = np.empty([DATA_IN_MEMORY_TIME, n, 6+2])
    positionCarbon = \
        np.empty([DATA_IN_MEMORY_TIME, int(n/(x_range.shape[0] - 1)), 6+2])
    positionHelium = np.empty([DATA_IN_MEMORY_TIME, n, 6+2])
    # Сразу крупные частицы
    num = 0  # номер частицы
    num = 0  # номер углерода
    carbon_num = 0  # номер углерода
    time = 0  # позиция по временной шкале с учетом цикличности заполнения
    # TODO оптимизировать получения случвела Максвелла
    electron_randomizer = MWranomizer(n=1, m=ELECTRONS_MASS, T=TEMPERATURE)
    carbon_randomizer = MWranomizer(n=1, m=CARBONS_MASS, T=TEMPERATURE)
    helium_randomizer = MWranomizer(n=1, m=HELIUMS_MASS, T=HELIUMS_TEMPERATURE)
    csdu, esdu, hsdu = CARBON_SPEED_DIMENSIONLESS_UNIT, \
        ELECTRON_SPEED_DIMENSIONLESS_UNIT, HELLIUM_SPEED_DIMENSIONLESS_UNIT
    big_carbon_radius = \
        (hi * CARBONS_RADIUS ** 3 * CARBONS_NUMBER) ** (1.0 / 3.0)
    big_carbon_charge = CARBONS_CHARGE * CARBONS_NUMBER
    big_electron_radius = \
        (hi * ELECTRONS_RADIUS ** 3 * ELECTRONS_NUMBER)**(1.0 / 3.0)
    big_electron_charge = ELECTRONS_CHARGE * ELECTRONS_NUMBER
    big_helium_radius = \
        (hi * HELIUMS_RADIUS ** 3 * HELIUMS_NUMBER) ** (1.0 / 3.0)
    big_helium_charge = HELIUMS_CHARGE * HELIUMS_NUMBER
    for y, j in zip(y_range[:-1], range(y_range.shape[0] - 1)):
        for z, k in zip(z_range[:-1], range(z_range.shape[0] - 1)):
            # Распределение углерода
            cell = (0, y, z)
            x_big, y_big, z_big = \
                spread_position(cell, CARBONS_NUMBER)
            x_speed, y_speed, z_speed = \
                spread_speed(carbon_randomizer, dimensionless=csdu)
            for v, l in zip([0, y_big, z_big, big_carbon_radius, big_carbon_charge, x_speed, y_speed, z_speed], range(positionCarbon.shape[2])):
                positionCarbon[time][carbon_num][l] = v
            carbon_num += 1
            for x, i in zip(x_range[:-1], range(x_range.shape[0] - 1)):
                cell = (x, y, z)
                print(cell)
                # Распределение электронов
                x_big, y_big, z_big = \
                    spread_position(cell, ELECTRONS_NUMBER)
                x_speed, y_speed, z_speed = \
                    spread_speed(electron_randomizer, dimensionless=esdu)
                for v, l in zip([x_big, y_big, z_big, big_electron_radius, big_electron_charge, x_speed, y_speed, z_speed], range(positionElectron.shape[2])):
                    positionElectron[time][num][l] = v
                # Распределение гелия
                x_big, y_big, z_big = \
                    spread_position(cell, HELIUMS_NUMBER)
                x_speed, y_speed, z_speed = \
                    spread_speed(helium_randomizer, dimensionless=hsdu)
                for v, l in zip([x_big, y_big, z_big, big_helium_radius, big_helium_charge, x_speed, y_speed, z_speed], range(positionHelium.shape[2])):
                    positionHelium[time][num][l] = v
                num += 1

    # Speed distribution 
    # begin_speed_distribution_data = [Mnuc*np.sqrt(positionCarbon[time][i][5]**2 + positionCarbon[time][i][6]**2 + positionCarbon[time][i][7]**2)  for i in range(positionCarbon.shape[1])]
    begin_speed_distribution_data = [positionCarbon[time][i][5]  for i in range(positionCarbon.shape[1])]
    end_speed_distribution_data = []
    
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
        while (time < MODELING_TIME):
            curr_time = p_time(time)
            # Заряд в узлах
            electron_charge_grid = np.zeros([x_range.shape[0], y_range.shape[0], z_range.shape[0]])
            carbon_charge_grid = np.zeros([x_range.shape[0], y_range.shape[0], z_range.shape[0]])
            helium_charge_grid = np.zeros([x_range.shape[0], y_range.shape[0], z_range.shape[0]])
            positions = [positionElectron, positionCarbon, positionHelium]
            grids = [electron_charge_grid, carbon_charge_grid, helium_charge_grid]
            names = ['positionElectron', 'positionCarbon', 'positionHelium']
            for grid, position, name in zip(grids, positions, names):
                for num in range(position.shape[1]):
                    # print(name)
                    x_big, y_big, z_big, _, charge, _, _, _ = position[curr_time][num]
                    i, j, k = \
                        int(x_big/X_STEP), int(y_big/Y_STEP), int(z_big/Z_STEP)
                    x, y, z = \
                        i*X_STEP, j*Y_STEP, k*Z_STEP
                    patch = spread_charge((x, y, z), (x_big, y_big, z_big), charge)
                    for p in patch:
                        grid[i+p[0]][j+p[1]][k+p[2]] += p[3]
            # Граничные условия для потенциала
            n = (x_range.shape[0]+2, y_range.shape[0]+2, z_range.shape[0]+2)
            phi = POTENTIAL_BOUND_VALUE
            ecg, ccg, hcg = \
                electron_charge_grid, carbon_charge_grid, helium_charge_grid
            prev_phi, next_phi, ro = \
                make_boundary_conditions(phi, n, ecg, ccg, hcg)
            # Метод установления
            prev_phi, next_phi = \
                potential_establish_method_cuda(prev_phi, next_phi, ro, epsilon=ESTABLISHING_METHOD_ACCURACY)

            # Расчет напряженности
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

            # Расчет напряженности действующей на частицу
            electron_tension = np.empty([positionElectron.shape[1], 3])
            carbon_tension = np.empty([positionCarbon.shape[1], 3])
            helium_tension = np.empty([positionHelium.shape[1], 3])
            particles = [positionElectron, positionCarbon, positionHelium]
            tensions = [electron_tension, carbon_tension, helium_tension]
            for position, p in zip(particles, range(len(tensions))):
                n = position.shape
                for num in range(n[1]):
                    x_big, y_big, z_big = \
                        get_component(position[curr_time][num])
                    i, j, k = \
                        int(x_big/X_STEP), int(y_big/Y_STEP), int(z_big/Z_STEP)
                    x, y, z = i*X_STEP, j*Y_STEP, k*Z_STEP
                    tension = spread_tension((x, y, z), (x_big, y_big, z_big), intensity)
                    for t in range(3):
                        tensions[p][num][t] = tension[t]
            # Решение дифуров
            t = np.linspace(0, TIME_STEP, 2)
            curr_time = p_time(time)
            for num in range(positionCarbon.shape[1]):
                x_big, y_big, z_big = \
                    get_component(positionCarbon[curr_time][num])
                i, j, k = \
                    int(x_big/X_STEP), int(y_big/Y_STEP), int(z_big/Z_STEP)
                # speeds = getSpeedProjection(positionCarbon[curr_time][num])
                speeds = get_component(positionCarbon[curr_time][num], b=5)

                for l in range(positionCarbon.shape[2]):
                    positionCarbon[p_next_time(time)][num][l] = positionCarbon[curr_time][num][l]
                for dim in range(3):
                    v = speeds[dim]
                    r = positionCarbon[curr_time][num][dim]/Ml
                    E = carbon_tension[num][dim]/Mee
                    y0 = [r, v, E]
                    res = odeint(f, y0, t)
                    # if num == 0:
                    #     print(res)
                    positionCarbon[p_next_time(time)][num][dim] = res[-1][0]*Ml
                    positionCarbon[p_next_time(time)][num][5+dim] = res[-1][1]
            for num in range(positionElectron.shape[1]):
                for l in range(positionElectron.shape[2]):
                    positionElectron[p_next_time(time)][num][l] = positionElectron[curr_time][num][l]
            for num in range(positionHelium.shape[1]):
                for l in range(positionHelium.shape[2]):
                    positionHelium[p_next_time(time)][num][l] = positionHelium[curr_time][num][l]

            time += 1
            print('time = {} '.format(time))

            make_tracks_plot_file(prefix, positionCarbon, time)
            make_intent_plot_file(prefix, intent, time)
            make_potential_plot_file(prefix, next_phi, time)
            make_intensity_plot_file(prefix, intensity, time)
            make_tension_plot_file(prefix, carbon_tension, time)

            curr_time = p_prev_time(time)
            for p, l in zip(listen_particles, range(len(listen_particles))):
                plot_data[l] += [positionCarbon[curr_time][p].copy()]
    except IndexError:
        print('IndexError')
        pass
    except KeyboardInterrupt:
        print('KeyboardInterrupt')
        pass

    # directory = make_dir(prefix, 'graphs')

    # speed distributionsudo numbersElectron            
    end_speed_distribution_data = [Mnuc*np.sqrt(positionCarbon[curr_time][i][5]**2 + positionCarbon[curr_time][i][6]**2 + positionCarbon[curr_time][i][7]**2)  for i in range(positionCarbon.shape[1])]
    print (begin_speed_distribution_data)
    print (end_speed_distribution_data)
    make_distribution_plot_file(prefix, begin_speed_distribution_data, end_speed_distribution_data, time)
    make_speed_position_plot_file(prefix, plot_data, listen_particles, time)


def make_dir_prefix():
    prefix = './picts/picts_cuda_{}'.format(time.strftime('%Y%m%d%H%M%S'))
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
    main(prefix)
