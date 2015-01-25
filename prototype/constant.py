# -*- coding: utf-8 -*-

__author__ = 'fedul0x'

import numpy as np
from numpy.random import rand


# Такты моделирования
# modelingTime = 100
# modeling_time = 1200
# modeling_time = 85
modeling_time = 1200
in_memory_time = 5
# Электростатическая постоянная
epsilon0 = 8.8541878E-12
# Число частиц в ячейке
# numbersElectron = 1E+10
# numbersCarbon = 1E+16
# numbersHelium = 1E+10
numbersElectron = 100
numbersCarbon = 100
numbersHelium = 100

# Массы частиц
massElectron = 9.11E-31
massCarbon = 12.011 * 1.67E-27 - 9.11E-31
massHelium = 4.002 * 1.67E-27 - 9.11E-31
# Заряды частиц
chargeElectron = -1.602176565E-19
chargeCarbon = 1.602176565E-19
chargeHelium = 1.602176565E-19
# Радиусы частиц
radiusElectron = 2.81794E-15
radiusCarbon = 31E-12
radiusHelium = 91E-12
# deltaT = 1.0E-12
# deltaT = 1.0E-2
deltaT = 1.0E-8
# deltaT = 1
# Время метода установления
# deltaT_ = 1.0E-12
deltaT_ = 1.75E-10

# Скорость выгорания анода
burningSpeedAnode = 0.44E-5

temperature = 4200
temperatureHelium = 293
amperage = 150

kb = 1.38E-23
Boltzmann_constant = 1.380648813E-23

maxSpeedElectron = np.sqrt(2 * kb * temperature / massElectron)
deltaSpeedElectron = 0.1*maxSpeedElectron
maxSpeedHelium = np.sqrt(2 * kb * temperatureHelium / massHelium)
deltaSpeedHelium = 0.1*maxSpeedHelium
maxSpeedCarbon = np.sqrt(2 * kb * temperature / massCarbon)
deltaSpeedCarbon = 0.1*maxSpeedCarbon

# Размеры сетки
xDimensionGrid = 1.0E-3
yDimensionGrid = 1.0E-2
zDimensionGrid = 1.0E-2
# число шагов по сетке для крупных частиц
xNumberStepGrid = 5  #im
yNumberStepGrid = 5  #km
zNumberStepGrid = 5  #lm
xInitNumberStepGrid = 50  #imm
yInitNumberStepGrid = 50  #kmm
zInitNumberStepGrid = 50  #lmm
# шаг по сетке
xStepGrid = xDimensionGrid / xNumberStepGrid  #hx
yStepGrid = yDimensionGrid / yNumberStepGrid  #hy
zStepGrid = zDimensionGrid / zNumberStepGrid  #hz
# шаг по сетке(для начального положения частиц)
xInitStepGrid = xDimensionGrid / xInitNumberStepGrid  #hxx
yInitStepGrid = yDimensionGrid / yInitNumberStepGrid  #hyy
zInitStepGrid = zDimensionGrid / zInitNumberStepGrid  #hzz

# Коэффициенты скорости
# Kvu2 = 0.005
# Kvu1 = (3-2*Kvu2**2)
# TODO Исправить коэффициенты
Kvu2 = 0.47619047619047616
Kvu1 = 0.04761904761904767

# Коэффициенты скорости для перехода к безразмерным величинам
Mnue = np.sqrt(2*kb*temperature/massElectron)
Mnuc = np.sqrt(2*kb*temperature/massCarbon)
Mnuh = np.sqrt(2*kb*temperatureHelium/massHelium)

# Коэффициенты для перехода к безразмерному пространству
Ml = np.sqrt(kb*epsilon0*temperature)/np.abs(chargeElectron)/np.sqrt(numbersCarbon)

# Коэффициенты для перехода к безразмерной напряженности
Mee = 8090.769400

# Метод установления
# epsilon = 0.001
epsilon = 0.01

# Mve = 3.5671E+5

# ЧТо такое hi
hi = 10



