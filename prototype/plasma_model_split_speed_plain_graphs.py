# -*- coding: utf-8 -*-

import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# import constant
from constant import *

# def dif(arr, i, j, k, index, step):
#     val = 0
#     if index == 0:
#         val = (arr[i+1][j][k] - 2*arr[i][j][k] + arr[i-1][j][k])/step**2
#     if index == 1:
#         val = (arr[i][j+1][k] - 2*arr[i][j][k] + arr[i][j-1][k])/step**2
#     if index == 2:
#         val = (arr[i][j][k+1] - 2*arr[i][j][k] + arr[i][j][k-1])/step**2
#     return val

def plain_plot(carbon_plt):
    colors = ['red', 'green', 'blue']
    for g, c in zip(carbon_plt, range(len(carbon_plt))):
        plt.plot([i[0] for i in g], [i[1] for i in g], color=colors[c])
        plt.title("y from x c = {}".format(c))
        # plt.show();
        plt.savefig("y from x c={}.png".format(c))
        plt.clf()

    for g, c in zip(carbon_plt, range(len(carbon_plt))):
        plt.plot([i[0] for i in g], [i[2] for i in g], color=colors[c])
        plt.title("z from x c = {}".format(c))
        # plt.show();
        plt.savefig("z from x c={}.png".format(c))
        plt.clf()

    for g, c in zip(carbon_plt, range(len(carbon_plt))):
        plt.plot([i[1] for i in g], [i[2] for i in g], color=colors[c])
        break
        plt.title("z from y c = {}".format(c))
        # plt.show();
        plt.savefig("z from y c={}.png".format(c))
        plt.clf()

    for g, c in zip(carbon_plt, range(len(carbon_plt))):
        plt.plot([i[0] for i in g], [np.sqrt(i[5]**2+i[6]**2+i[7]**2)*np.sqrt(3)*Mnuc for i in g], color=colors[c])
    plt.title("v from x")
    plt.savefig("v from x.png")
    # plt.show();
    plt.clf()

    for g, c in zip(carbon_plt, range(len(carbon_plt))):
        plt.plot([i[1] for i in g], [np.sqrt(i[5]**2+i[6]**2+i[7]**2)*np.sqrt(3)*Mnuc for i in g], color=colors[c])
    plt.title("v from y")
    plt.savefig("v from y.png")
    # plt.show();
    plt.clf()


    for g, c in zip(carbon_plt, range(len(carbon_plt))):
        plt.plot([i[2] for i in g], [np.sqrt(i[5]**2+i[6]**2+i[7]**2)*np.sqrt(3)*Mnuc for i in g], color=colors[c])
    plt.title("v from z")
    plt.savefig("v from z.png")
    # plt.show();
    plt.clf()


def getSpeedProjection(particle):
    # speed = particle[5]
    # devider = np.sqrt(3)
    # return (speed*Kvu1/devider, speed*Kvu2/devider, speed*Kvu2/devider)
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


    # chargeGridElectron[i][j][k] += charge * oppositeX * oppositeY * oppositeZ / v
    # chargeGridElectron[i + 1][j][k] += charge * neighborX * oppositeY * oppositeZ / v
    # chargeGridElectron[i][j + 1][k] += charge * oppositeX * neighborY * oppositeZ / v
    # chargeGridElectron[i][j][k + 1] += charge * oppositeX * oppositeY * neighborZ / v
    # chargeGridElectron[i + 1][j + 1][k] += charge * neighborX * neighborY * oppositeZ / v
    # chargeGridElectron[i][j + 1][k + 1] += charge * oppositeX * neighborY * neighborZ / v
    # chargeGridElectron[i + 1][j][k + 1] += charge * neighborX * oppositeY * neighborZ / v
    # chargeGridElectron[i + 1][j + 1][k + 1] += charge * neighborX * neighborY * neighborZ / v
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


def main():
    # unsigned int r = 2000000
    # srand(r)

    corpusculNumbers = 1000
    corpusculNumbersPerCell = corpusculNumbers / (xNumberStepGrid * yNumberStepGrid * zNumberStepGrid)
    corpusculNumbersDelta = 0.1 * corpusculNumbersPerCell
    k = 0
    kkk = []
    # print(corpusculNumbersPerCell)
    xRange = np.linspace(0, xDimensionGrid, xNumberStepGrid + 1)
    yRange = np.linspace(0, yDimensionGrid, yNumberStepGrid + 1)
    zRange = np.linspace(0, zDimensionGrid, zNumberStepGrid + 1)
    # print(xRange, xStepGrid)


    # ЧТо такое hi
    hi = 10
    # positionElectron = np.empty([xRange.shape[0] - 1, yRange.shape[0] - 1, zRange.shape[0] - 1, 3])
    # positionCarbon = np.empty([xRange.shape[0] - 1, yRange.shape[0] - 1, zRange.shape[0] - 1, 3])
    # positionHelium = np.empty([xRange.shape[0] - 1, yRange.shape[0] - 1, zRange.shape[0] - 1, 3])
    num = (xRange.shape[0] - 1) * (yRange.shape[0] - 1) * (zRange.shape[0] - 1) 
    # time, i, x, y, z, radius, charge, speedx, , speedy, speedz 
    positionElectron = np.empty([modelingTime, num, 6+2])
    positionCarbon = np.empty([modelingTime, num, 6+2])
    positionHelium = np.empty([modelingTime, num, 6+2])
    # заряды в узлах сетки
    chargeBigElectron = np.empty([xRange.shape[0] - 1, yRange.shape[0] - 1, zRange.shape[0] - 1])
    chargeBigCarbon = np.empty([xRange.shape[0] - 1, yRange.shape[0] - 1, zRange.shape[0] - 1])
    chargeBigHelium = np.empty([xRange.shape[0] - 1, yRange.shape[0] - 1, zRange.shape[0] - 1])

    # сhargeGridElectron = np.array(xRange, yRange, zRange)
    chargeGridElectron = np.zeros([xRange.shape[0], yRange.shape[0], zRange.shape[0]])
    chargeGridCarbon = np.zeros([xRange.shape[0], yRange.shape[0], zRange.shape[0]])
    chargeGridHelium = np.zeros([xRange.shape[0], yRange.shape[0], zRange.shape[0]])
    # Сразу крупные частицы
    num = 0
    time = 0
    for x, i in zip(xRange[:-1], range(xRange.shape[0] - 1)):
        for y, j in zip(yRange[:-1], range(yRange.shape[0] - 1)):
            for z, k in zip(zRange[:-1], range(zRange.shape[0] - 1)):
                # Распределение электронов
                cn = int(np.round(corpusculNumbersPerCell + (rand() * 2 * corpusculNumbersDelta - corpusculNumbersDelta)))
                xBig, yBig, zBig, speed = spreadParticle((x, y, z), (xStepGrid, yStepGrid, zStepGrid), (maxSpeedElectron, deltaSpeedElectron), cn)
                speed = speed/Mnue
                speed = 0.0
                radiusBigElectron = (hi * radiusElectron ** 3 * cn) ** (1.0 / 3.0)
                # print(radiusBigElectron)
                # return
                chargeBigElectron = chargeElectron * cn
                for v, l in zip([xBig, yBig, zBig, radiusBigElectron, chargeBigElectron, speed, speed, speed], range(positionElectron.shape[2])):
                    positionElectron[time][num][l] = v
                patch = spreadCharge((x, y, z), (xStepGrid, yStepGrid, zStepGrid), (xBig, yBig, zBig), chargeBigElectron)
                for p in patch:
                    chargeGridElectron[i+p[0]][j+p[1]][k+p[2]] += p[3]
                # Распределение углерода
                cn = int(
                    np.round(corpusculNumbersPerCell + (rand() * 2 * corpusculNumbersDelta - corpusculNumbersDelta)))
                xBig, yBig, zBig, speed = spreadParticle((0, y, z), (xStepGrid, yStepGrid, zStepGrid), (maxSpeedCarbon, deltaSpeedCarbon), cn)
                speed = speed/Mnuc
                speed = 0.0
                radiusBigCarbon = (hi * radiusCarbon ** 3 * cn) ** (1.0 / 3.0) 
                chargeBigCarbon = chargeCarbon * cn
                for v, l in zip([0, yBig, zBig, radiusBigCarbon, chargeBigCarbon, speed, speed, speed], range(positionCarbon.shape[2])):
                    positionCarbon[time][num][l] = v

                patch = spreadCharge((x, y, z), (xStepGrid, yStepGrid, zStepGrid), (xBig, yBig, zBig), chargeBigCarbon)
                for p in patch:
                    chargeGridCarbon[i+p[0]][j+p[1]][k+p[2]] += p[3]

                # Распределение гелия
                cn = int(
                    np.round(corpusculNumbersPerCell + (rand() * 2 * corpusculNumbersDelta - corpusculNumbersDelta)))
                xBig, yBig, zBig, speed = spreadParticle((x, y, z), (xStepGrid, yStepGrid, zStepGrid), (maxSpeedHelium, deltaSpeedHelium), cn)
                speed = speed/Mnuh
                speed = 0.0
                radiusBigHelium = (hi * radiusHelium ** 3 * cn) ** (1.0 / 3.0)
                chargeBigHelium = chargeHelium * cn
                # positionHelium[num][0], positionHelium[num][1], positionHelium[num][2], positionHelium[num][3], positionHelium[num][4], positionHelium[num][5] = xBig, yBig, zBig, radiusBigHelium, chargeBigHelium, speed
                for v, l in zip([xBig, yBig, zBig, radiusBigHelium, chargeBigHelium, speed, speed, speed], range(positionHelium.shape[2])):
                    positionHelium[time][num][l] = v

                patch = spreadCharge((x, y, z), (xStepGrid, yStepGrid, zStepGrid), (xBig, yBig, zBig), chargeBigHelium)
                for p in patch:
                    chargeGridHelium[i+p[0]][j+p[1]][k+p[2]] += p[3]

                num += 1

    # Граничные условия для потенциала
    phiphi = 50
    phiphi_ = 49.9999995
    # phiphi_ = 0
    n = (chargeGridElectron.shape[0]+2, chargeGridElectron.shape[1]+2, chargeGridElectron.shape[2]+2)
    prevPhi = np.zeros([n[0], n[1], n[2]])
    nextPhi = np.zeros([n[0], n[1], n[2]])
    # intr = np.linspace(0, 50, n[0])
    # for j in range(n[1]):
    #     for k in range(n[2]):
    #         for i in range(n[0]):
    #                 prevPhi[i][j][k] = intr[-i-1] #+ np.random.random()*2
    #                 nextPhi[i][j][k] = intr[-i-1] #+ np.random.random()*2
    for i in range(n[0]):
        for j in range(n[1]):
            for k in range(n[2]):
                if (i * j * k == 0) or (i == n[0] - 1) or (j == n[1] - 1) or (k == n[2] - 1):
                # if (i == 0) or (i == n[0] - 1):
                    prevPhi[i][j][k] = phiphi  #50
                    nextPhi[i][j][k] = phiphi  #50
                else:
                    prevPhi[i][j][k] = phiphi_  #phisl
                    nextPhi[i][j][k] = phiphi_  #phisl
    #             # if (i == n[0] - 1):
    #             #     prevPhi[i][j][k] = 50.5#1.3*phiphi  #50
    #             #     nextPhi[i][j][k] = 50.5#1.3*phiphi #50

    # MODELING CYCLE BEGIN
    num = 0
    time = 0
    allTime = 0
    crashesElectron = set()
    lastCrashesElectron = set()
    crashesHelium = set()
    lastCrashesHelium = set()
    carbon_plt = [[], [], []]
    # print(positionElectron[0])
    while (time < modelingTime):
        # Заряд в узлах
        chargeGridElectron = np.zeros([xRange.shape[0], yRange.shape[0], zRange.shape[0]])
        chargeGridCarbon = np.zeros([xRange.shape[0], yRange.shape[0], zRange.shape[0]])
        chargeGridHelium = np.zeros([xRange.shape[0], yRange.shape[0], zRange.shape[0]])
        # n = (xRange.shape[0] - 1) * (yRange.shape[0] - 1) * (zRange.shape[0] - 1)
        for num in range(positionElectron.shape[1]):
            xBig, yBig, zBig, _, charge, _, _, _ = positionElectron[time][num]
            i, j, k = int(xBig/xStepGrid), int(yBig/yStepGrid), int(zBig/zStepGrid)
            if i>=chargeGridElectron.shape[0]-1 or j>=chargeGridElectron.shape[1]-1 or k>=chargeGridElectron.shape[2]-1 : 
                plain_plot(carbon_plt)
                return
            x, y, z = i*xStepGrid, j*yStepGrid, k*zStepGrid
            patch = spreadCharge((x, y, z), (xStepGrid, yStepGrid, zStepGrid), (xBig, yBig, zBig), charge)
            for p in patch:
                chargeGridElectron[i+p[0]][j+p[1]][k+p[2]] += p[3]
        for num in range(positionCarbon.shape[1]):
            xBig, yBig, zBig, _, charge, _, _, _ = positionCarbon[time][num]
            i, j, k = int(xBig/xStepGrid), int(yBig/yStepGrid), int(zBig/zStepGrid)
            if i>=chargeGridElectron.shape[0]-1 or j>=chargeGridElectron.shape[1]-1 or k>=chargeGridElectron.shape[2]-1 : 
                plain_plot(carbon_plt)
                return
            x, y, z = i*xStepGrid, j*yStepGrid, k*zStepGrid
            patch = spreadCharge((x, y, z), (xStepGrid, yStepGrid, zStepGrid), (xBig, yBig, zBig), charge)
            for p in patch:
                chargeGridCarbon[i+p[0]][j+p[1]][k+p[2]] += p[3]
        for num in range(positionHelium.shape[1]):
            xBig, yBig, zBig, _, charge, _, _, _ = positionHelium[time][num]
            i, j, k = int(xBig/xStepGrid), int(yBig/yStepGrid), int(zBig/zStepGrid)
            if i>=chargeGridElectron.shape[0]-1 or j>=chargeGridElectron.shape[1]-1 or k>=chargeGridElectron.shape[2]-1 : 
                plain_plot(carbon_plt)
                return
            x, y, z = i*xStepGrid, j*yStepGrid, k*zStepGrid
            patch = spreadCharge((x, y, z), (xStepGrid, yStepGrid, zStepGrid), (xBig, yBig, zBig), charge)
            for p in patch:
                chargeGridHelium[i+p[0]][j+p[1]][k+p[2]] += p[3]
        # print(chargeGridElectron)
        # print(chargeGridCarbon)
        # print(chargeGridHelium)
        # for x, i in zip(xRange[:-1], range(xRange.shape[0] - 1)):
        #     for y, j in zip(yRange[:-1], range(yRange.shape[0] - 1)):
        #         for z, k in zip(zRange[:-1], range(zRange.shape[0] - 1)):
        # print(prevPhi)

        # Метод установления
        n = (chargeGridElectron.shape[0]+2, chargeGridElectron.shape[1]+2, chargeGridElectron.shape[2]+2)
        prevSum, nextSum = epsilon, 3 * epsilon
        while np.abs(prevSum - nextSum) > epsilon:
            prevSum, nextSum = 0, 0
            for i in range(1, n[0] - 1):
                for j in range(1, n[1] - 1):
                    for k in range(1, n[2] - 1):
                        ro = chargeGridElectron[i-1][j-1][k-1] + chargeGridCarbon[i-1][j-1][k-1] + chargeGridHelium[i-1][j-1][k-1]
                        # print('\n {} {} {}'.format(i, j, k))
                        # print('ro = {}'.format(ro))
                        # print('dif1 = {}'.format(dif(prevPhi[i + 1][j][k], prevPhi[i][j][k], prevPhi[i - 1][j][k], xStepGrid)))
                        # print('dif2 = {}'.format(dif(prevPhi[i][j + 1][k], prevPhi[i][j][k], prevPhi[i][j - 1][k], yStepGrid)))
                        # print('dif3 = {}'.format(dif(prevPhi[i][j][k + 1], prevPhi[i][j][k], prevPhi[i][j][k - 1], zStepGrid)))
                        # print('dif1+dif2+dif3 + 4piRO = {}'.format((
                        #     dif(prevPhi[i + 1][j][k], prevPhi[i][j][k], prevPhi[i - 1][j][k], xStepGrid) +
                        #     dif(prevPhi[i][j + 1][k], prevPhi[i][j][k], prevPhi[i][j - 1][k], yStepGrid) +
                        #     dif(prevPhi[i][j][k + 1], prevPhi[i][j][k], prevPhi[i][j][k - 1], zStepGrid) +
                        #     4 * np.pi * ro)))
                        # print("rhs = {}".format(deltaT_ * (
                        #     dif(prevPhi[i + 1][j][k], prevPhi[i][j][k], prevPhi[i - 1][j][k], xStepGrid) +
                        #     dif(prevPhi[i][j + 1][k], prevPhi[i][j][k], prevPhi[i][j - 1][k], yStepGrid) +
                        #     dif(prevPhi[i][j][k + 1], prevPhi[i][j][k], prevPhi[i][j][k - 1], zStepGrid) +
                        #     4 * np.pi * ro)))
                        # input()
                        nextPhi[i][j][k] = deltaT_ * (
                            dif(prevPhi[i + 1][j][k], prevPhi[i][j][k], prevPhi[i - 1][j][k], xStepGrid) +
                            dif(prevPhi[i][j + 1][k], prevPhi[i][j][k], prevPhi[i][j - 1][k], yStepGrid) +
                            dif(prevPhi[i][j][k + 1], prevPhi[i][j][k], prevPhi[i][j][k - 1], zStepGrid) +
                            4 * np.pi * ro) + prevPhi[i][j][k]
                        prevSum += prevPhi[i][j][k]
                        nextSum += nextPhi[i][j][k]
            for i in range(chargeGridElectron.shape[0] ):
                for j in range(chargeGridElectron.shape[1] ):
                    for k in range(chargeGridElectron.shape[2] ):
                        prevPhi[i][j][k] = nextPhi[i][j][k]
            print('{}>{}'.format(np.abs(prevSum - nextSum), epsilon))
        print(nextPhi)

        # Расчет напряженности
        n = nextPhi.shape
        intensity = np.zeros([n[0]-2, n[1]-2, n[2]-2, 3])
        for i in range(1, n[0]-1):
            for j in range(1, n[1]-1):
                for k in range(1, n[2]-1):
                    intensity[i-1][j-1][k-1][0] = (nextPhi[i - 1][j][k] - nextPhi[i + 1][j][k]) / 2 / xStepGrid
                    intensity[i-1][j-1][k-1][1] = (nextPhi[i][j - 1][k] - nextPhi[i][j + 1][k]) / 2 / yStepGrid
                    intensity[i-1][j-1][k-1][2] = (nextPhi[i][j][k - 1] - nextPhi[i][j][k + 1]) / 2 / zStepGrid

        # Расчет напряженности действующей на частицу
        n = positionElectron.shape
        tensionCorpusculActingElectron = np.empty([xRange.shape[0] - 1, yRange.shape[0] - 1, zRange.shape[0] - 1, 3])
        for num in range(positionElectron.shape[1]):
            # print(time, num)
            xBig, yBig, zBig = positionElectron[time][num][0], positionElectron[time][num][1], positionElectron[time][num][2]
            i, j, k = int(xBig/xStepGrid), int(yBig/yStepGrid), int(zBig/zStepGrid)
            x, y, z = i*xStepGrid, j*yStepGrid, k*zStepGrid
            # print(xBig, yBig, zBig)
            tension = spreadTension((x, y, z), (xStepGrid, yStepGrid, zStepGrid), (xBig, yBig, zBig), intensity)
            for t in range(3):
                tensionCorpusculActingElectron[i][j][k][t] = tension[t]

        n = positionCarbon.shape
        tensionCorpusculActingCarbon = np.empty([xRange.shape[0] - 1, yRange.shape[0] - 1, zRange.shape[0] - 1, 3])
        for num in range(positionCarbon.shape[1]):
            xBig, yBig, zBig = positionCarbon[time][num][0], positionCarbon[time][num][1], positionCarbon[time][num][2]
            i, j, k = int(xBig/xStepGrid), int(yBig/yStepGrid), int(zBig/zStepGrid)
            x, y, z = i*xStepGrid, j*yStepGrid, k*zStepGrid
            # print(xBig, yBig, zBig)
            tension = spreadTension((x, y, z), (xStepGrid, yStepGrid, zStepGrid), (xBig, yBig, zBig), intensity)
            for t in range(3):
                tensionCorpusculActingCarbon[i][j][k][t] = tension[t]

        n = positionHelium.shape
        tensionCorpusculActingHelium = np.empty([xRange.shape[0] - 1, yRange.shape[0] - 1, zRange.shape[0] - 1, 3])
        for num in range(positionHelium.shape[1]):
            xBig, yBig, zBig = positionHelium[time][num][0], positionHelium[time][num][1], positionHelium[time][num][2]
            i, j, k = int(xBig/xStepGrid), int(yBig/yStepGrid), int(zBig/zStepGrid)
            x, y, z = i*xStepGrid, j*yStepGrid, k*zStepGrid
            # print(xBig, yBig, zBig)
            tension = spreadTension((x, y, z), (xStepGrid, yStepGrid, zStepGrid), (xBig, yBig, zBig), intensity)
            for t in range(3):
                tensionCorpusculActingHelium[i][j][k][t] = tension[t]

        # time += 1

        # Стокновения углерода и электронов
        # pe = positionElectron
        # pc = positionCarbon
        # lastCrashesElectron = crashesElectron
        # crashesElectron = set()
        # for e in range(pe.shape[1]):
        #     for c in range(pc.shape[1]):
        #         # print("{} <= {}".format(np.sqrt((pe[time][e][0] - pc[time][c][0] )**2 + (pe[time][e][1] - pc[time][c][1] )**2 + (pe[time][e][2] - pc[time][c][2] )**2),  pe[time][e][3] + pc[time][c][3]))
        #         # print("{} > {}".format(np.sqrt((pe[time][e][0] - pc[time][c][0] )**2 + (pe[time][e][1] - pc[time][c][1] )**2 + (pe[time][e][2] - pc[time][c][2] )**2), np.abs(pe[time][e][3] - pc[time][c][3])))
        #         if np.sqrt((pe[time][e][0] - pc[time][c][0] )**2 + (pe[time][e][1] - pc[time][c][1] )**2 + (pe[time][e][2] - pc[time][c][2] )**2)  <= pe[time][e][3] + pc[time][c][3]: # and \
        #             # np.sqrt((pe[time][e][0] - pc[time][c][0] )**2 + (pe[time][e][1] - pc[time][c][1] )**2 + (pe[time][e][2] - pc[time][c][2] )**2) > np.abs(pe[time][e][3] - pc[time][c][3]):

        #             spd1 = [(2 * massElectron * pe[time][e][5]  + pc[time][c][5]* (massCarbon - massElectron))/(massCarbon + massElectron) for p in range(5, 8)]
        #             spd2 = [(2 * massCarbon * pc[time][c][5]  + pe[time][e][5]* (massElectron - massCarbon))/(massCarbon + massElectron) for p in range(5, 8)]

        #             pc[time][c][5], pc[time][c][6], pc[time][c][7] = spd1[0], spd1[1], spd1[2]
        #             pe[time][e][5], pe[time][e][6], pe[time][e][7] = spd2[0], spd2[1], spd2[2]
        #             crashesElectron.add(e)
        # ph = positionHelium
        # pc = positionCarbon
        # lastCrashesHelium = crashesHelium
        # crashesHelium = set()
        # for h in range(ph.shape[1]):
        #     for c in range(pc.shape[1]):
        #         # print("{} <= {}".format(np.sqrt((ph[time][h][0] - pc[time][c][0] )**2 + (ph[time][h][1] - pc[time][c][1] )**2 + (ph[time][h][2] - pc[time][c][2] )**2),  ph[time][h][3] + pc[time][c][3]))
        #         # print("{} > {}".format(np.sqrt((ph[time][h][0] - pc[time][c][0] )**2 + (ph[time][h][1] - pc[time][c][1] )**2 + (ph[time][h][2] - pc[time][c][2] )**2), np.abs(ph[time][h][3] - pc[time][c][3])))
        #         if np.sqrt((ph[time][h][0] - pc[time][c][0] )**2 + (ph[time][h][1] - pc[time][c][1] )**2 + (ph[time][h][2] - pc[time][c][2] )**2)  <= ph[time][h][3] + pc[time][c][3]: #and \
        #             # np.sqrt((ph[time][h][0] - pc[time][c][0] )**2 + (ph[time][h][1] - pc[time][c][1] )**2 + (ph[time][h][2] - pc[time][c][2] )**2) > np.abs(ph[time][h][3] - pc[time][c][3]):

        #             spd1 = [(2 * massHelium * ph[time][h][5]  + pc[time][c][5]* (massCarbon - massHelium))/(massCarbon + massHelium) for p in range(5, 8)]
        #             spd2 = [(2 * massCarbon * pc[time][c][5]  + ph[time][h][5]* (massHelium - massCarbon))/(massCarbon + massHelium) for p in range(5, 8)]

        #             pc[time][c][5], pc[time][c][6], pc[time][c][7] = spd1[0], spd1[1], spd1[2]
        #             ph[time][h][5], ph[time][h][6], ph[time][h][7] = spd2[0], spd2[1], spd2[2]
        #             crashesHelium.add(h)

        from scipy.integrate import odeint
        def f(y, t):
            delta = 1.0
            Z = 1.0
            # epsilon = 8.854187817E-12
            epsilonC = 1.0
            r,v = y[0], y[1]
            f0 = np.sqrt(delta)*v
            # f1 = np.sqrt(delta)*Z/2/epsilon0*E_
            f1 = np.sqrt(delta)*Z/2/epsilonC*E_
            return [f0, f1]

        t = np.linspace(0, 2*deltaT, 20)
        # t = np.linspace(0, 2*deltaT*1000000*1000000, 20)
        # pceusoTime = (time+1)
        pceusoTime = (time+1) % (modelingTime-1)
        for num in range(positionCarbon.shape[1]):
            xBig, yBig, zBig, _, _, _, _, _ = positionCarbon[time][num]
            i, j, k = int(xBig/xStepGrid), int(yBig/yStepGrid), int(zBig/zStepGrid)
            speeds = getSpeedProjection(positionCarbon[time][num])

            for l in range(positionCarbon.shape[2]):
                positionCarbon[pceusoTime][num][l] = positionCarbon[time][num][l]
            positionCarbon[pceusoTime][num][5] = 0.0
            positionCarbon[pceusoTime][num][6] = 0.0
            positionCarbon[pceusoTime][num][7] = 0.0
            for dim in range(3):
                v = speeds[dim]
                r = positionCarbon[time][num][dim]/Ml
                E_ = tensionCorpusculActingCarbon[i][j][k][dim]/Mee
                y0 = [r, v]
                res = odeint(f, y0, t)
                positionCarbon[pceusoTime][num][dim] = res[-1][0]*Ml
                positionCarbon[pceusoTime][num][5+dim] = res[-1][1]
                if num == 0:
                    print("coord {}: time={} y0={} E= {} ^^^^".format(dim, time, y0, E_))
                # print(res)
            # positionCarbon[pceusoTime][num][5] = np.sqrt(positionCarbon[pceusoTime][num][5])
        for num in range(positionElectron.shape[1]):
            for l in range(positionElectron.shape[2]):
                positionElectron[pceusoTime][num][l] = positionElectron[time][num][l]
        for num in range(positionHelium.shape[1]):
            for l in range(positionHelium.shape[2]):
                positionHelium[pceusoTime][num][l] = positionHelium[time][num][l]
        # print("Ml = {}".format(Ml))


        # time += 1
        time = pceusoTime


        # for plotting
        allTime += 1


        # from mpl_toolkits.mplot3d import proj3d
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # ax.set_xlim3d([0, xDimensionGrid])
        # ax.set_ylim3d([0, yDimensionGrid])
        # ax.set_zlim3d([0, zDimensionGrid])
        # tensionCorpusculActing = tensionCorpusculActingCarbon
        # p = []
        # import matplotlib.patches as mpatches
        # labels = [[[i for i in range(positionElectron.shape[0])] for j in range(positionElectron.shape[0])] for k in range(positionElectron.shape[0])]
        # for num in range(positionCarbon.shape[0]):
        #     # xBig, yBig, zBig = positionElectron[time][num][0], positionElectron[time][num][1], positionElectron[time][num][2]
        #     xBig, yBig, zBig = positionCarbon[time][num][0], positionCarbon[time][num][1], positionCarbon[time][num][2]
        #     # xBig, yBig, zBig = positionHelium[time][num][0], positionHelium[time][num][1], positionHelium[time][num][2]
        #     i, j, k = int(xBig/xStepGrid), int(yBig/yStepGrid), int(zBig/zStepGrid)
        #     x, y, z = i*xStepGrid, j*yStepGrid, k*zStepGrid
        #     # print(xBig, yBig, zBig)
        #     ax.scatter(xBig, yBig, zBig, color='green')  #, s=int( radiusBig*1E+14)+35)
        #     # if (i%2 == j%2 == k%2 == 0):
        #     #     x2, y2, _ = proj3d.proj_transform(xBig, yBig, zBig, ax.get_proj())
        #     #     str_label = '{:5.2f},{:5.2f},{:5.2f}'.format(tensionCorpusculActing[i][j][k][0],tensionCorpusculActing[i][j][k][1],tensionCorpusculActing[i][j][k][2])
        #     #     str_point = '{},{},{}'.format(i, j, k)
        #     #     p = p + [mpatches.Patch(color='green', label='({}) - ({})'.format(str_point, str_label))]
        #     #     print(str_point)
        #         # labels[i][j][k] = plt.annotate (str_point, xy=(x2, y2), xytext=(-7, 7),
        #         #     textcoords='offset points', ha='right', va='bottom',
        #         #     bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
        #         #     arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

        # # def update_position(e):
        # #     for num in range(positionCarbon.shape[0]):
        # #         xBig, yBig, zBig = positionElectron[num][0], positionElectron[num][1], positionElectron[num][2]
        # #         i, j, k = int(xBig/xStepGrid), int(yBig/yStepGrid), int(zBig/zStepGrid)
        # #         x, y, z = i*xStepGrid, j*yStepGrid, k*zStepGrid
        # #         if (not(i%2 == j%2 == k%2 == 0)):
        # #             continue
        # #         x2, y2, _ = proj3d.proj_transform(xBig, yBig, zBig, ax.get_proj())
        # #         labels[i][j][k].xy = x2,y2
        # #         labels[i][j][k].update_positions(fig.canvas.renderer)

        # # fig.canvas.draw()
        # # fig.canvas.mpl_connect('button_release_event', update_position)
        # plt.title('time = {}'.format(allTime))
        # # plt.legend(handles=p)
        # # plt.show()
        # plt.savefig("./picts/carbon_3d_{0:04d}.png".format(allTime))
        # plt.clf()
        





        # from mpl_toolkits.mplot3d import proj3d
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # ax.set_xlim3d([0, xDimensionGrid])
        # ax.set_ylim3d([0, yDimensionGrid])
        # ax.set_zlim3d([0, zDimensionGrid])
        # tensionCorpusculActing = tensionCorpusculActingCarbon
        # p = []
        # import matplotlib.patches as mpatches
        # labels = [[[i for i in range(positionElectron.shape[0])] for j in range(positionElectron.shape[0])] for k in range(positionElectron.shape[0])]
        # flags = [[False for ii in range(6)] for jj in range(6)]
        # for num in range(positionCarbon.shape[0]):
        #     # xBig, yBig, zBig = positionElectron[time][num][0], positionElectron[time][num][1], positionElectron[time][num][2]
        #     xBig, yBig, zBig = positionCarbon[time][num][0], positionCarbon[time][num][1], positionCarbon[time][num][2]
        #     # xBig, yBig, zBig = positionHelium[time][num][0], positionHelium[time][num][1], positionHelium[time][num][2]
        #     i, j, k = int(xBig/xStepGrid), int(yBig/yStepGrid), int(zBig/zStepGrid)
        #     x, y, z = i*xStepGrid, j*yStepGrid, k*zStepGrid
        #     ax.scatter(xBig, yBig, zBig, color='red')  #, s=int( radiusBig*1E+14)+35)
        #     if (j == k and not flags[j][k]):
        #         flags[j][k] = True
        #         x2, y2, _ = proj3d.proj_transform(xBig, yBig, zBig, ax.get_proj())
        #         str_label = '{:5.2f},{:5.2f},{:5.2f}'.format(positionCarbon[time][num][5], positionCarbon[time][num][6], positionCarbon[time][num][7])
        #         # str_label = '{:5.2f},{:5.2f},{:5.2f}'.format(tensionCorpusculActing[i][j][k][0],tensionCorpusculActing[i][j][k][1],tensionCorpusculActing[i][j][k][2])
        #         # str_point = '{},{},{}'.format(i, j, k)
        #         # p = p + [mpatches.Patch(color='green', label='({}) - ({})'.format(str_point, str_label))]
        #         # print(str_point)
        #         labels[i][j][k] = plt.annotate (str_label, xy=(x2, y2), xytext=(-15, 15),
        #             textcoords='offset points', ha='right', va='bottom',
        #             bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
        #             arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

        # for num, cnt in zip(crashesElectron | lastCrashesElectron, range(len(crashesElectron | lastCrashesElectron))):
        #     xBig, yBig, zBig = positionElectron[time][num][0], positionElectron[time][num][1], positionElectron[time][num][2]
        #     ax.scatter(xBig, yBig, zBig, color='green')  #, s=int( radiusBig*1E+14)+35)
        #     x2, y2, _ = proj3d.proj_transform(xBig, yBig, zBig, ax.get_proj())
        #     str_label = '{:5.2f},{:5.2f},{:5.2f}'.format(positionElectron[time][num][5], positionElectron[time][num][6], positionElectron[time][num][7])
        #     plt.annotate (str_label, xy=(x2, y2), xytext=(+100+20*cnt, -50),
        #         textcoords='offset points', ha='right', va='bottom',
        #         bbox=dict(boxstyle='round,pad=0.5', fc='green', alpha=0.5),
        #         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

        # for num in crashesHelium | lastCrashesHelium:
        #     xBig, yBig, zBig = positionHelium[time][num][0], positionHelium[time][num][1], positionHelium[time][num][2]
        #     ax.scatter(xBig, yBig, zBig, color='blue')  #, s=int( radiusBig*1E+14)+35)
        #     x2, y2, _ = proj3d.proj_transform(xBig, yBig, zBig, ax.get_proj())
        #     str_label = '{:5.2f},{:5.2f},{:5.2f}'.format(positionHelium[time][num][5], positionHelium[time][num][6], positionHelium[time][num][7])
        #     plt.annotate (str_label, xy=(x2, y2), xytext=(+150, 50),
        #         textcoords='offset points', ha='right', va='bottom',
        #         bbox=dict(boxstyle='round,pad=0.5', fc='blue', alpha=0.5),
        #         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
        # print(crashesElectron | lastCrashesElectron)
        # print(crashesElectron, lastCrashesElectron)
        # print(crashesHelium | lastCrashesHelium)
        # print(crashesHelium, lastCrashesHelium)

        # # fig.canvas.draw()
        # # fig.canvas.mpl_connect('button_release_event', update_position)
        # plt.title('time = {} \n phix = 50.0..0.0, electron crashes: {} '.format(allTime, len(crashesElectron)))
        # # plt.legend(handles=p)
        # plt.show()
        # # plt.savefig("./picts_crashes/carbon_3d_speed_{0:04d}.png".format(allTime))
        # plt.clf()




        # plt.xlim(0, xDimensionGrid)
        # plt.ylim(0, yDimensionGrid)
        # tensionCorpusculActing = tensionCorpusculActingCarbon
        # p = []
        # import matplotlib.patches as mpatches
        # labels = [[[i for i in range(positionElectron.shape[0])] for j in range(positionElectron.shape[0])] for k in range(positionElectron.shape[0])]
        # for num in range(positionCarbon.shape[1]):
        #     # xBig, yBig, zBig = positionElectron[time][num][0], positionElectron[time][num][1], positionElectron[time][num][2]
        #     xBig, yBig, zBig = positionCarbon[time][num][0], positionCarbon[time][num][1], positionCarbon[time][num][2]
        #     # xBig, yBig, zBig = positionHelium[time][num][0], positionHelium[time][num][1], positionHelium[time][num][2]
        #     i, j, k = int(xBig/xStepGrid), int(yBig/yStepGrid), int(zBig/zStepGrid)
        #     x, y, z = i*xStepGrid, j*yStepGrid, k*zStepGrid
        #     # print(xBig, yBig, zBig)
        #     plt.scatter(xBig, yBig, color='green')  #, s=int( radiusBig*1E+14)+35)
        #     pc = positionCarbon
        #     plt.plot([pc[t][num][0] for t in range(pceusoTime+1)], [pc[t][num][1] for t in range(pceusoTime+1)])

        
        # plt.title('time = {}'.format(allTime))
        # # plt.legend(handles=p)

        # # plt.show()

        # plt.savefig("./picts/carbon_2d_{0:04d}.png".format(allTime))
        # plt.clf()

        carbon_plt[0] += [positionCarbon[time][0]]
        carbon_plt[1] += [positionCarbon[time][int(positionCarbon.shape[1]/2)]]
        carbon_plt[2] += [positionCarbon[time][positionCarbon.shape[1]-1]]




if __name__ == '__main__':
    main()