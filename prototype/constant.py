# -*- coding: utf-8 -*-
import numpy as np

__author__ = 'fedul0x'

# Параметры моделирования
FAST_ESTABLISHING_METHOD = True
# Время моделирования
MODELING_TIME = 1200
DATA_IN_MEMORY_TIME = 3
# Электростатическая постоянная
# epsilon0 = 8.8541878E-12
ELECTROSTATIC_CONSTANT = 8.8541878E-12
# Постоянная Больцмана
# kb = 1.38E-23
BOLTZMANN_CONSTANT = 1.380648813E-23
# Число Авогадро
AVOGADRO_CONSTANT = 6.0221412927E+23
# Число частиц в ячейке
# numbersElectron = 1E+10
# numbersCarbon = 1E+16
# numbersHelium = 1E+10
ELECTRONS_NUMBER = 1E+12
CARBONS_NUMBER = 1E+12
HELIUMS_NUMBER = 1E+6
# Массы частиц
ELECTRONS_MASS = 9.11E-31
CARBONS_MASS = 12.011 * 1.67E-27 - 9.11E-31
HELIUMS_MASS = 4.002 * 1.67E-27 - 9.11E-31
# Заряды частиц
ELECTRONS_CHARGE = -1.602176565E-19
CARBONS_CHARGE = 1.602176565E-19
HELIUMS_CHARGE = 1.602176565E-19
# Радиусы частиц
ELECTRONS_RADIUS = 2.81794E-15
CARBONS_RADIUS = 91E-12
HELIUMS_RADIUS = 31E-12
# CARBONS_RADIUS, HELIUMS_RADIUS = HELIUMS_RADIUS, CARBONS_RADIUS
# deltaT = 1.0E-12
# deltaT = 1.0E-2
# deltaT = 1.0E-8
TIME_STEP = 1.0E-2
# TIME_STEP = 1.0E-1
# deltaT = 1.0E-8
# deltaT = 1
# Скорость выгорания анода
ANODE_BURNING_SPEED = 0.44E-5
# Температура процесса
# TEMPERATURE = 4200
TEMPERATURE = 10000
HELIUMS_TEMPERATURE = 293
AMPERAGE = 150
# Параметры распределения скорости по Максвеллу
ELECTRONS_MAX_SPEED = np.sqrt(2*BOLTZMANN_CONSTANT*TEMPERATURE/ELECTRONS_MASS)
ELECTRONS_DELTA_SPEED = 0.1*ELECTRONS_MAX_SPEED
HELIUMS_MAX_SPEED = \
    np.sqrt(2*BOLTZMANN_CONSTANT*HELIUMS_TEMPERATURE/HELIUMS_MASS)
HELIUMS_DELTA_SPEED = 0.1*HELIUMS_MAX_SPEED
CARBONS_MAX_SPEED = np.sqrt(2*BOLTZMANN_CONSTANT*TEMPERATURE/CARBONS_MASS)
CARBONS_DELTA_SPEED = 0.1*CARBONS_MAX_SPEED
# Размеры сетки
X_DIMENSION_GRID = 1.0E-3
Y_DIMENSION_GRID = 1.0E-2
Z_DIMENSION_GRID = 1.0E-2
# Число шагов по сетке для крупных частиц
X_STEP_NUMBER_GRID = 50
Y_STEP_NUMBER_GRID = 50
Z_STEP_NUMBER_GRID = 50
# X_STEP_NUMBER_GRID = 10
# Y_STEP_NUMBER_GRID = 10
# Z_STEP_NUMBER_GRID = 10
# Шаг по сетке
X_STEP = X_DIMENSION_GRID / X_STEP_NUMBER_GRID
Y_STEP = Y_DIMENSION_GRID / Y_STEP_NUMBER_GRID
Z_STEP = Z_DIMENSION_GRID / Z_STEP_NUMBER_GRID
# Начальное значение времени метода установления
FAKE_TIME_STEP = 1.75E-9
while (FAKE_TIME_STEP/X_STEP**2 + FAKE_TIME_STEP/Y_STEP**2 + FAKE_TIME_STEP/Z_STEP**2 > 0.5):
	FAKE_TIME_STEP /= 2
# print('FAKE_TIME_STEP = {}'.format(FAKE_TIME_STEP))
# Начальное значение потенциала
POTENTIAL_BOUND_VALUE = 50
# Коэффициенты скорости для перехода к безразмерным величинам
ELECTRON_SPEED_DIMENSIONLESS_UNIT = \
    np.sqrt(2*BOLTZMANN_CONSTANT*TEMPERATURE/ELECTRONS_MASS)
CARBON_SPEED_DIMENSIONLESS_UNIT = \
    np.sqrt(2*BOLTZMANN_CONSTANT*TEMPERATURE/CARBONS_MASS)
HELLIUM_SPEED_DIMENSIONLESS_UNIT = \
    np.sqrt(2*BOLTZMANN_CONSTANT*HELIUMS_TEMPERATURE/HELIUMS_MASS)
Mnue = np.sqrt(2*BOLTZMANN_CONSTANT*TEMPERATURE/ELECTRONS_MASS)
Mnuc = np.sqrt(2*BOLTZMANN_CONSTANT*TEMPERATURE/CARBONS_MASS)
Mnuh = np.sqrt(2*BOLTZMANN_CONSTANT*HELIUMS_TEMPERATURE/HELIUMS_MASS)

# Коэффициенты для перехода к безразмерному пространству
SPACE_DIMENSIONLESS_UNIT = \
    np.sqrt(BOLTZMANN_CONSTANT*ELECTROSTATIC_CONSTANT*TEMPERATURE) \
    / np.abs(ELECTRONS_CHARGE) / np.sqrt(CARBONS_NUMBER)
Ml = np.sqrt(BOLTZMANN_CONSTANT*ELECTROSTATIC_CONSTANT*TEMPERATURE) \
    / np.abs(ELECTRONS_CHARGE) / np.sqrt(CARBONS_NUMBER)
# Коэффициенты для перехода к безразмерной напряженности
INTENSITY_DIMENSIONLESS_UNIT = 8090.769400
Mee = 8090.769400
# Метод установления
# epsilon = 0.001
# epsilon = 0.01
ESTABLISHING_METHOD_ACCURACY = 0.01
# Mve = 3.5671E+5
# ЧТо такое hi
hi = 0.5  #0..1



