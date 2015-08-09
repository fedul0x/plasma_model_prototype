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
from collision import *


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

class MWranomizer2:

    def __init__(self, n, m, T):
        ''' Maxwell idstribution randomizer '''
        self.a = 0
        self.n, self.m, self.T = n, m, T

    def rvs(self, size=1):
        return scipy.stats.maxwell.rvs(scale=CARBONS_MAX_SPEED/1.5957691216057308, size=size)

def get_component(particle, b=0, n=3):
    return [particle[i] for i in range(b, b+n)]


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

def spread_speed(randomizer, dimensionless=1):
    speed = randomizer.rvs(size=1)[0]/dimensionless
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
    electron = np.empty([DATA_IN_MEMORY_TIME, n, 6+2])
    carbon = \
        np.empty([DATA_IN_MEMORY_TIME, int(n/(x_range.shape[0] - 1)), 6+2])
    helium = np.empty([DATA_IN_MEMORY_TIME, n, 6+2])
    # Сразу крупные частицы
    num = 0  # номер частицы
    num = 0  # номер углерода
    carbon_num = 0  # номер углерода
    time = 0  # позиция по временной шкале с учетом цикличности заполнения
    # TODO оптимизировать получения случвела Максвелла
    print('Particle distribution')
    electron_randomizer = MWranomizer2(n=1, m=ELECTRONS_MASS, T=TEMPERATURE)
    carbon_randomizer = MWranomizer2(n=1, m=CARBONS_MASS, T=TEMPERATURE)
    helium_randomizer = MWranomizer2(n=1, m=HELIUMS_MASS, T=HELIUMS_TEMPERATURE)
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
            for v, l in zip([0, y_big, z_big, big_carbon_radius, big_carbon_charge, x_speed, y_speed, z_speed], range(carbon.shape[2])):
                carbon[time][carbon_num][l] = v
            carbon_num += 1
            for x, i in zip(x_range[:-1], range(x_range.shape[0] - 1)):
                cell = (x, y, z)
                # print(j, k, i)
                # Распределение электронов
                x_big, y_big, z_big = \
                    spread_position(cell, ELECTRONS_NUMBER)
                x_speed, y_speed, z_speed = \
                    spread_speed(electron_randomizer, dimensionless=esdu)
                for v, l in zip([x_big, y_big, z_big, big_electron_radius, big_electron_charge, x_speed, y_speed, z_speed], range(electron.shape[2])):
                    electron[time][num][l] = v
                # Распределение гелия
                x_big, y_big, z_big = \
                    spread_position(cell, HELIUMS_NUMBER)
                x_speed, y_speed, z_speed = \
                    spread_speed(helium_randomizer, dimensionless=hsdu)
                for v, l in zip([x_big, y_big, z_big, big_helium_radius, big_helium_charge, x_speed, y_speed, z_speed], range(helium.shape[2])):
                    helium[time][num][l] = v
                num += 1

    # Speed distribution 
    begin_speed_distribution_data = [Mnuc*np.sqrt(carbon[time][i][5]**2 + carbon[time][i][6]**2 + carbon[time][i][7]**2)  for i in range(carbon.shape[1])]
    # begin_speed_distribution_data = [carbon[time][i][5]  for i in range(carbon.shape[1])]
    end_speed_distribution_data = []
    
    # MODELING CYCLE BEGIN
    num = 0 # номер частицы
    time = 0 # абсолютная позиция по временной шкале
    crashesElectron = []
    lastCrashesElectron = []
    crashesHelium = []
    lastCrashesHelium = []
    prev_phi, next_phi = None, None
    # Для итоговых графиков
    listen_particles = [0, int(carbon.shape[1]/2), carbon.shape[1]-1]
    plot_data = [[[]] for _ in listen_particles]
    typeI, typeII, typeIII = [], [], []
    max_distances = []
    try:
        while (time < MODELING_TIME):
            curr_time = p_time(time)
            # Заряд в узлах
            print('Calc charge in cells')
            electron_charge_grid = np.zeros([x_range.shape[0], y_range.shape[0], z_range.shape[0]])
            carbon_charge_grid = np.zeros([x_range.shape[0], y_range.shape[0], z_range.shape[0]])
            helium_charge_grid = np.zeros([x_range.shape[0], y_range.shape[0], z_range.shape[0]])
            positions = [electron, carbon, helium]
            grids = [electron_charge_grid, carbon_charge_grid, helium_charge_grid]
            names = ['electron', 'carbon', 'helium']
            for grid, position, name in zip(grids, positions, names):
                n = grid.shape
                for num in range(position.shape[1]):
                    # print(name)
                    x_big, y_big, z_big, _, charge, _, _, _ = position[curr_time][num]
                    try:
                        i, j, k = \
                            int(x_big/X_STEP), int(y_big/Y_STEP), int(z_big/Z_STEP)
                    except Exception:
                        print('{}/{}, {}/{}, {}/{}'.format(x_big, X_STEP, y_big, Y_STEP, z_big, Z_STEP))
                    i, j, k = \
                        int(x_big/X_STEP), int(y_big/Y_STEP), int(z_big/Z_STEP)

                    if name == 'carbon':
                        if (i<0) or (j<0) or (k<0) or (i>n[0]-2) or (j>n[1]-2) or (k>n[2]-2):
                            cell = (0, np.random.choice(y_range[:-1]), np.random.choice(z_range[:-1]))
                            x_big, y_big, z_big = \
                                spread_position(cell, CARBONS_NUMBER)
                            x_speed, y_speed, z_speed = \
                                spread_speed(carbon_randomizer, dimensionless=csdu)
                            for v, l in zip([0, y_big, z_big, big_carbon_radius, big_carbon_charge, x_speed, y_speed, z_speed], range(carbon.shape[2])):
                                position[curr_time][num][l] = v
                            x_big, y_big, z_big, _, charge, _, _, _ = position[curr_time][num]
                            i, j, k = \
                                int(x_big/X_STEP), int(y_big/Y_STEP), int(z_big/Z_STEP)
                            if num in listen_particles:
                                plot_data[listen_particles.index(num)] += []

                            # print(i, j, k)

                        # print('int(x_big/X_STEP)={}/{}, int(y_big/Y_STEP)=={}/{}, int(z_big/Z_STEP)=={}/{}'.format(x_big, X_STEP, y_big, Y_STEP, z_big, Z_STEP))
                    x, y, z = \
                        i*X_STEP, j*Y_STEP, k*Z_STEP
                    patch = spread_charge((x, y, z), (x_big, y_big, z_big), charge)
                    # print(name)
                    for p in patch:
                        grid[i+p[0]][j+p[1]][k+p[2]] += p[3]
            # Граничные условия для потенциала
            n = (x_range.shape[0]+2, y_range.shape[0]+2, z_range.shape[0]+2)
            phi = POTENTIAL_BOUND_VALUE
            ecg, ccg, hcg = \
                electron_charge_grid, carbon_charge_grid, helium_charge_grid
            # TODO easying
            if FAST_ESTABLISHING_METHOD:
                if None == prev_phi and None == next_phi:
                    prev_phi, next_phi, ro = \
                        make_boundary_conditions(phi, n, ecg, ccg, hcg)
                else:
                    print('ECONOMIC')
                    _, _, ro = \
                        make_boundary_conditions(phi, n, ecg, ccg, hcg)
            else:
                prev_phi, next_phi, ro = \
                    make_boundary_conditions(phi, n, ecg, ccg, hcg)
            # Метод установления
            print('Establishing method')
            prev_phi, next_phi = \
                potential_establish_method_cuda(prev_phi, next_phi, ro, epsilon=ESTABLISHING_METHOD_ACCURACY)
                # potential_establish_method(prev_phi, next_phi, ro, epsilon=ESTABLISHING_METHOD_ACCURACY)

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
            electron_tension = np.empty([electron.shape[1], 3])
            carbon_tension = np.empty([carbon.shape[1], 3])
            helium_tension = np.empty([helium.shape[1], 3])
            particles = [electron, carbon, helium]
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

            print('Crashing')
            crashesCarbon, crashesElectron, crashesHelium = [], [], []
            if time != 0:
                res = find_collision(time, electron, carbon, helium, max_distances)
                crashesCarbon, crashesElectron, crashesHelium = res[0]
                typeI.extend(res[1][0])
                typeII.extend(res[1][1])
                typeIII.extend(res[1][2])

            # Решение дифуров
            t = np.linspace(0, TIME_STEP, 2)
            curr_time = p_time(time)
            for num in range(carbon.shape[1]):
                x_big, y_big, z_big = \
                    get_component(carbon[curr_time][num])
                i, j, k = \
                    int(x_big/X_STEP), int(y_big/Y_STEP), int(z_big/Z_STEP)
                # speeds = getSpeedProjection(carbon[curr_time][num])
                speeds = get_component(carbon[curr_time][num], b=5)

                for l in range(carbon.shape[2]):
                    carbon[p_next_time(time)][num][l] = carbon[curr_time][num][l]
                for dim in range(3):
                    # import pdb
                    # pdb.set_trace()
                    v = speeds[dim]
                    r = carbon[curr_time][num][dim]/SPACE_DIMENSIONLESS_UNIT
                    E = carbon_tension[num][dim]/INTENSITY_DIMENSIONLESS_UNIT
                    y0 = [r, v, E]
                    res = odeint(f, y0, t)
                    # if num == 0:
                    #     print(res)
                    carbon[p_next_time(time)][num][dim] = res[-1][0]*SPACE_DIMENSIONLESS_UNIT
                    carbon[p_next_time(time)][num][5+dim] = res[-1][1]
                # print(carbon[p_next_time(time)][num][0], carbon[p_next_time(time)][num][1], carbon[p_next_time(time)][num][2])
                # print(carbon[p_next_time(time)][num][5], carbon[p_next_time(time)][num][6], carbon[p_next_time(time)][num][7])
            for num in range(electron.shape[1]):
                for l in range(electron.shape[2]):
                    electron[p_next_time(time)][num][l] = electron[curr_time][num][l]
            for num in range(helium.shape[1]):
                for l in range(helium.shape[2]):
                    helium[p_next_time(time)][num][l] = helium[curr_time][num][l]


            # # Стокновения углерода и электронов
            # pe = electron
            # ph = helium
            # pc = carbon
            # # lastCrashesElectron = crashesElectron
            # crashesElectron = []
            # crashesHelium = []
            # print('1')

            # for c in range(pc.shape[1]):
            #     for e in range(pe.shape[1]):
            #         print(e, c)
            #         if np.sqrt((pe[curr_time][e][0] - pc[curr_time][c][0] )**2 + (pe[curr_time][e][1] - pc[curr_time][c][1] )**2 + (pe[curr_time][e][2] - pc[curr_time][c][2] )**2)  <= pe[curr_time][e][3] + pc[curr_time][c][3]  and \
            #             np.sqrt((pe[curr_time][e][0] - pc[curr_time][c][0] )**2 + (pe[curr_time][e][1] - pc[curr_time][c][1] )**2 + (pe[curr_time][e][2] - pc[curr_time][c][2] )**2) > np.abs(pe[curr_time][e][3] - pc[curr_time][c][3]):

            #             crashesElectron += [pe[curr_time][e].copy()]
            #             # spd1 = [(2 * massElectron * pe[curr_time][e][p]  + pc[curr_time][c][p]* (massCarbon - massElectron))/(massCarbon + massElectron) for p in range(5, 8)]
            #             # spd2 = [(2 * massCarbon * pc[curr_time][c][p]  + pe[curr_time][e][p]* (massElectron - massCarbon))/(massCarbon + massElectron) for p in range(5, 8)]

            #             # pc[curr_time][c][5], pc[curr_time][c][6], pc[curr_time][c][7] = spd1[0], spd1[1], spd1[2]
            #             # pe[curr_time][e][5], pe[curr_time][e][6], pe[curr_time][e][7] = spd2[0], spd2[1], spd2[2]
            #     for h in range(ph.shape[1]):
            #         if np.sqrt((ph[curr_time][h][0] - pc[curr_time][c][0] )**2 + (ph[curr_time][h][1] - pc[curr_time][c][1] )**2 + (ph[curr_time][h][2] - pc[curr_time][c][2] )**2)  <= ph[curr_time][h][3] + pc[curr_time][c][3] and \
            #             np.sqrt((ph[curr_time][h][0] - pc[curr_time][c][0] )**2 + (ph[curr_time][h][1] - pc[curr_time][c][1] )**2 + (ph[curr_time][h][2] - pc[curr_time][c][2] )**2) > np.abs(ph[curr_time][h][3] - pc[curr_time][c][3]):

            #             crashesHelium += [ph[curr_time][h][:]]
            #             # spd1 = [(2 * massHelium * ph[curr_time][h][p]  + pc[curr_time][c][p]* (massCarbon - massHelium))/(massCarbon + massHelium) for p in range(5, 8)]
            #             # spd2 = [(2 * massCarbon * pc[curr_time][c][p]  + ph[curr_time][h][p]* (massHelium - massCarbon))/(massCarbon + massHelium) for p in range(5, 8)]
            #             # pc[curr_time][c][5], pc[curr_time][c][6], pc[curr_time][c][7] = spd1[0], spd1[1], spd1[2]
            #             # ph[curr_time][h][5], ph[curr_time][h][6], ph[curr_time][h][7] = spd2[0], spd2[1], spd2[2]
            # # pc = carbon
            # # print('1.5')
            # # # lastCrashesHelium = crashesHelium
            # #     for c in range(pc.shape[1]):
            # #         print(h, c)
            # # print('2')
            # Новые столкновения

            time += 1
            print('time = {} '.format(time))

            # directory = make_dir(prefix, 'graphs')
            # plt.plot([i[0] for i in max_distances], color='red', label='x')
            # plt.plot([0, len(max_distances)], [X_STEP_COLLISION, X_STEP_COLLISION], linestyle=':', lw=6, color='red')
            # plt.plot([i[1] for i in max_distances], color='green', label='y')
            # plt.plot([0, len(max_distances)], [Y_STEP_COLLISION, Y_STEP_COLLISION], linestyle=':', lw=6, color='green')
            # plt.plot([i[2] for i in max_distances], color='blue', label='z')
            # plt.plot([0, len(max_distances)], [Z_STEP_COLLISION, Z_STEP_COLLISION], linestyle=':', lw=6, color='blue')
            # plt.legend()
            # plt.savefig("{}/max_distances".format(directory))
            # plt.clf()
            # plt.close()
            
            
            # make_tracks_plot_file(prefix, carbon, time)
            make_tracks_plot_file_with_crashes(prefix, carbon, set(crashesElectron), set(crashesHelium), set(crashesCarbon), time)
            # make_inten_plot_file(prefix, inten, time)
            # make_potential_plot_file(prefix, next_phi, time)
            # make_intensity_plot_file(prefix, intensity, time)
            # make_tension_plot_file(prefix, carbon_tension, time)

            curr_time = p_prev_time(time)
            for p, l in zip(listen_particles, range(len(listen_particles))):
                plot_data[l][-1] += [carbon[curr_time][p].copy()]
            msg = 'LENGTH {}'.format(len(plot_data))
            for q in plot_data:
                msg += '\n   {}'.format(len(q))
                for w in q:
                    msg += '\n     {}'.format(len(w))
                    # for e in w:
                    #     msg += '\n       {}'.format(len(e))
            print(msg)

            plt.clf()
            plt.close()
    # except IndexError:
    #     print('IndexError')
    #     pass
    except KeyboardInterrupt:
        print('KeyboardInterrupt')
        pass
    except OSError:
        print('OSError')
        pass

    # directory = make_dir(prefix, 'graphs')

    # speed distributionsudo numbersElectron            
    end_speed_distribution_data = [Mnuc*np.sqrt(carbon[curr_time][i][5]**2 + carbon[curr_time][i][6]**2 + carbon[curr_time][i][7]**2)  for i in range(carbon.shape[1])]
    print (begin_speed_distribution_data)
    print (end_speed_distribution_data)
    make_distribution_plot_file(prefix, begin_speed_distribution_data, end_speed_distribution_data)
    print('typeI: {} \ntypeII: {} \ntypeIII {}'.format(typeI, typeII, typeIII))
    make_collision_distribution_plot_file(prefix, typeI, typeII, typeIII)
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

def find_collision(time, electron, carbon, helium, max_distances):
    curr_time = p_time(time)
    prev_time = p_prev_time(time)

    max_distance_x, max_distance_y, max_distance_z = 0.0, 0.0, 0.0
    for i in range(carbon.shape[1]):
        value = np.abs(carbon[curr_time][i][0] - carbon[prev_time][i][0])
        if value > max_distance_x:
            max_distance_x = value
        value = np.abs(carbon[curr_time][i][1] - carbon[prev_time][i][1])
        if value > max_distance_y:
            max_distance_y = value
        value = np.abs(carbon[curr_time][i][2] - carbon[prev_time][i][2])
        if  value > max_distance_z:
            max_distance_z = value

    max_distances += [(max_distance_x, max_distance_y, max_distance_z)]
    # (max_distance_x, max_distance_y, max_distance_z)
    if max_distance_x == 0 or X_DIMENSION_GRID/max_distance_x > 100 or X_DIMENSION_GRID/max_distance_x < 60:
        X_DIMENSION_COLLISION = 100
    else:
        X_DIMENSION_COLLISION = int(X_DIMENSION_GRID/max_distance_x)

    if max_distance_y == 0 or Y_DIMENSION_GRID/max_distance_y > 100 or Y_DIMENSION_GRID/max_distance_y < 60:
        Y_DIMENSION_COLLISION = 100
    else:
        Y_DIMENSION_COLLISION = int(Y_DIMENSION_GRID/max_distance_y)
    
    if max_distance_z == 0 or Z_DIMENSION_GRID/max_distance_z > 100 or Z_DIMENSION_GRID/max_distance_z < 60:
        Z_DIMENSION_COLLISION = 100
    else:
        Z_DIMENSION_COLLISION = int(Z_DIMENSION_GRID/max_distance_z)
    
    # Z_DIMENSION_COLLISION = int(Y_DIMENSION_GRID/max_distance_y)
    # X_DIMENSION_COLLISION = int(Z_DIMENSION_GRID/max_distance_z)
    X_STEP_COLLISION = X_DIMENSION_GRID / X_DIMENSION_COLLISION
    Y_STEP_COLLISION = Y_DIMENSION_GRID / Y_DIMENSION_COLLISION
    Z_STEP_COLLISION = Z_DIMENSION_GRID / Z_DIMENSION_COLLISION
    msg = "MAX DIST BY COORDS: \n{}<{}, \n{}<{}, \n{}<{}"
    print(msg.format(max_distance_x, X_STEP_COLLISION, max_distance_y, Y_STEP_COLLISION, max_distance_z, Z_STEP_COLLISION))
    msg = "X_DIMENSION_COLLISION = {}, Y_DIMENSION_COLLISION = {}, Z_DIMENSION_COLLISION = {}"
    print(msg.format(X_DIMENSION_COLLISION, Y_DIMENSION_COLLISION, Z_DIMENSION_COLLISION))
    print('MAKE GRID')
    dims = (X_DIMENSION_COLLISION, Y_DIMENSION_COLLISION, Z_DIMENSION_COLLISION)
    steps = (X_STEP_COLLISION, Y_STEP_COLLISION, Z_STEP_COLLISION)

    electron_collision = make_collision_grid(X_DIMENSION_COLLISION, Y_DIMENSION_COLLISION, Z_DIMENSION_COLLISION)
    helium_collision = make_collision_grid(X_DIMENSION_COLLISION, Y_DIMENSION_COLLISION, Z_DIMENSION_COLLISION)
    carbon_collision = make_collision_grid(X_DIMENSION_COLLISION, Y_DIMENSION_COLLISION, Z_DIMENSION_COLLISION)
    electron_collision = distribute_by_grid(electron, electron_collision, curr_time, dims, steps)
    helium_collision = distribute_by_grid(helium, helium_collision, curr_time, dims, steps)
    carbon_collision = distribute_by_grid(carbon, carbon_collision, curr_time, dims, steps)
    crashesCarbon, crashesElectron, crashesHelium = [], [], []

    print('CRASH FINDING')
    count = 0
    typeI, typeII, typeIII = [], [], []
    for i in range(X_DIMENSION_COLLISION):
        for j in range(Y_DIMENSION_COLLISION):
            for k in range(Z_DIMENSION_COLLISION):
                for a in carbon_collision[i][j][k]:
                    xc, yc, zc, radiusc, _, vcx, vcy, vcz = carbon[curr_time][a] 
                    for b in electron_collision[i][j][k]:
                        xe, ye, ze, radiuse, _, vex, vey, vez = electron[curr_time][b]
                        # bp1 = get_component(carbon[p_time(time)][a])
                        # ep1 = get_component(carbon[p_next_time(time)][a])
                        # bp2 = get_component(carbon[p_time(time)][b])
                        # ep2 = get_component(carbon[p_next_time(time)][b])
                        # TODO remove 5*
                        # if np.sqrt((xc-xe)**2 + (yc-ye)**2 + (zc-ze)**2) < 5*(radiusc + radiuse):
                        bp1 = get_component(carbon[curr_time][a])
                        ep1 = get_component(carbon[prev_time][a])
                        bp2 = get_component(electron[curr_time][b])
                        ep2 = get_component(electron[prev_time][b])
                        vs = check_collision(bp1, ep1, bp2, ep2, radiusc, radiuse, TIME_STEP)
                        if vs[0]:
                            crashesElectron += [tuple(electron[curr_time][b][:])]
                            mc = CARBONS_MASS*radiusc**3 / hi / CARBONS_RADIUS**3
                            me = ELECTRONS_MASS*radiuse**3 / hi / ELECTRONS_RADIUS**3
                            electron[curr_time][b][5] = (2 * mc * vcx + vex*(me - mc))/(mc + me)
                            electron[curr_time][b][6] = (2 * mc * vcy + vey*(me - mc))/(mc + me)
                            electron[curr_time][b][7] = (2 * mc * vcz + vez*(me - mc))/(mc + me)
                            carbon[curr_time][a][5] = (2 * me * vex + vcx*(mc - me))/(mc + me)
                            carbon[curr_time][a][6] = (2 * me * vey + vcy*(mc - me))/(mc + me)
                            carbon[curr_time][a][7] = (2 * me * vez + vcz*(mc - me))/(mc + me)
                            # make_crash_track(prefix, bp1, ep1, bp2, ep2, radiusc, radiuse, vs, TIME_STEP, count)
                            count += 1
                            # print('E  ({:.4f}, {:.4f}, {:.4f})\n=> ({:.4f}, {:.4f}, {:.4f})'.format(vex, vey, vez, electron[curr_time][b][5], electron[curr_time][b][6], electron[curr_time][b][7]))
                            # print('C  ({:.4f}, {:.4f}, {:.4f})\n=> ({:.4f}, {:.4f}, {:.4f})'.format(vcx, vcy, vcz, carbon[curr_time][a][5], carbon[curr_time][a][6], carbon[curr_time][a][7]))
                    for b in helium_collision[i][j][k]:
                        xh, yh, zh, radiush, _, vhx, vhy, vhz = helium[curr_time][b]
                        # TODO remove 5*
                        # if np.sqrt((xc-xh)**2 + (yc-yh)**2 + (zc-zh)**2) < 5*(radiusc + radiush):
                        bp1 = get_component(carbon[curr_time][a])
                        ep1 = get_component(carbon[prev_time][a])
                        bp2 = get_component(helium[curr_time][b])
                        ep2 = get_component(helium[prev_time][b])
                        vs = check_collision(bp1, ep1, bp2, ep2, radiusc, radiush, TIME_STEP)
                        if vs[0]:
                            crashesHelium += [tuple(helium[curr_time][b][:])]
                            mc = CARBONS_MASS*radiusc**3 / hi / CARBONS_RADIUS**3
                            mh = HELIUMS_MASS*radiush**3 / hi / HELIUMS_RADIUS**3
                            helium[curr_time][b][5] = (2 * mc * vcx + vhx*(mh - mc))/(mc + mh)
                            helium[curr_time][b][6] = (2 * mc * vcy + vhy*(mh - mc))/(mc + mh)
                            helium[curr_time][b][7] = (2 * mc * vcz + vhz*(mh - mc))/(mc + mh)
                            carbon[curr_time][a][5] = (2 * mh * vhx + vcx*(mc - mh))/(mc + mh)
                            carbon[curr_time][a][6] = (2 * mh * vhy + vcy*(mc - mh))/(mc + mh)
                            carbon[curr_time][a][7] = (2 * mh * vhz + vcz*(mc - mh))/(mc + mh)
                            # make_crash_track(prefix, bp1, ep1, bp2, ep2, radiusc, radiush, vs, TIME_STEP, count)
                            count += 1
                            # print('H  ({:.4f}, {:.4f}, {:.4f})\n=> ({:.4f}, {:.4f}, {:.4f})'.format(vhx, vhy, vhz, helium[curr_time][b][5], helium[curr_time][b][6], helium[curr_time][b][7]))
                            # print('C  ({:.4f}, {:.4f}, {:.4f})\n=> ({:.4f}, {:.4f}, {:.4f})'.format(vcx, vcy, vcz, carbon[curr_time][a][5], carbon[curr_time][a][6], carbon[curr_time][a][7]))
                    for b in carbon_collision[i][j][k]:
                        xc_2, yc_2, zc_2, radiusc_2, _, vcx_2, vcy_2, vcz_2 = carbon[curr_time][b]
                        # TODO remove 45*
                        # if (a != b) and (np.sqrt((xc-xc_2)**2 + (yc-yc_2)**2 + (zc-zc_2)**2) < (radiusc + radiusc_2)):
                        bp1 = get_component(carbon[curr_time][a])
                        ep1 = get_component(carbon[prev_time][a])
                        bp2 = get_component(carbon[curr_time][b])
                        ep2 = get_component(carbon[prev_time][b])
                        vs = check_collision(bp1, ep1, bp2, ep2, radiusc, radiusc_2, TIME_STEP)
                        if (a != b) and vs[0]:
                            crashesCarbon += [tuple(carbon[curr_time][b][:])]
                            m_1 = CARBONS_MASS*radiusc**3 / hi / CARBONS_RADIUS**3
                            m_2 = CARBONS_MASS*radiusc_2**3 / hi / CARBONS_RADIUS**3
                            v_1 = CARBON_SPEED_DIMENSIONLESS_UNIT*np.sqrt(vcx**2 + vcy**2 + vcz**2)
                            v_2 = CARBON_SPEED_DIMENSIONLESS_UNIT*np.sqrt(vcx_2**2 + vcy_2**2 + vcz_2**2)
                            e_1 = m_1 * v_1**2 / 2 / CARBONS_NUMBER * AVOGADRO_CONSTANT
                            e_2 = m_2 * v_2**2 / 2 / CARBONS_NUMBER * AVOGADRO_CONSTANT
                            print('E1 = {} E2 = {} sum = {}'.format(e_1, e_2, e_1+e_2))
                            # 348000
                            if e_1 + e_2 > 348000:
                                typeI += [carbon[curr_time][b][0]]
                                continue
                                # 614000
                            if e_1 + e_2 > 614000:
                                typeII += [carbon[curr_time][b][0]]
                                continue
                                # 839000
                            if e_1 + e_2 > 839000:
                                typeIII += [carbon[curr_time][b][0]]
                                continue
                            # make_crash_tracks(prefix, bp1, ep1, bp2, ep2, radiusc, radiusc_2, vs, TIME_STEP, count)
                            count += 1
                            # e_1 = m_1 * v_1**2 / 2 * AVOGADRO_CONSTANT
                            # e_2 = m_2 * v_2**2 / 2 * AVOGADRO_CONSTANT
    return ((crashesCarbon, crashesElectron, crashesHelium), (typeI, typeII, typeIII))



if __name__ == '__main__':
    prefix = make_dir_prefix()
    main(prefix)