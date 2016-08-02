# -*- coding: utf-8 -*-
import numpy as np
import os.path

__author__ = 'fedul0x'

# Использовать ускоренный расчет эквипотенциальных поверхностей 
# (на основе предыдущих значений)
FAST_ESTABLISHING_METHOD = True

# Каталог с дампами, указать None при отсутвии оного
DUMP_FOLDER = '/home/fedul0x/tmp/dumps'
# DUMP_FOLDER = None

# Файл базы данных
DB_FILE = '/home/fedul0x/tmp/dumps/db.sqlite'

# Колчиство тактов моделирования
MODELING_TIME = 1200
# Длинна истории
DATA_IN_MEMORY_TIME = 3
# Количество отрывающихся слоев углерода
CARBON_LAYERS_NUMBER = 30

# Электростатическая постоянная (epsilon0)
ELECTROSTATIC_CONSTANT = 8.8541878E-12

# Постоянная Больцмана (kb)
BOLTZMANN_CONSTANT = 1.380648813E-23

# Число Авогадро
AVOGADRO_CONSTANT = 6.0221412927E+23

# Число частиц в ячейке
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

# Такт квантования
TIME_STEP = 1.0E-8

# Скорость выгорания анода
ANODE_BURNING_SPEED = 0.44E-5
# Температура процесса
TEMPERATURE = 4200
HELIUMS_TEMPERATURE = 293

# AMPERAGE = 150
AMPERAGE = 400

# Параметры распределения скорости по Максвеллу
ELECTRONS_MAX_SPEED = np.sqrt(2*BOLTZMANN_CONSTANT*TEMPERATURE/ELECTRONS_MASS)
ELECTRONS_DELTA_SPEED = 0.1*ELECTRONS_MAX_SPEED
HELIUMS_MAX_SPEED = \
    np.sqrt(2*BOLTZMANN_CONSTANT*HELIUMS_TEMPERATURE/HELIUMS_MASS)
HELIUMS_DELTA_SPEED = 0.1*HELIUMS_MAX_SPEED
CARBONS_MAX_SPEED = np.sqrt(2*BOLTZMANN_CONSTANT*TEMPERATURE/CARBONS_MASS)
CARBONS_DELTA_SPEED = 0.1*CARBONS_MAX_SPEED
# Дисперсия угла при распределении скорости по компонентам
SPEED_DISTRIBUTION_ANGLE_SIGMA = np.pi/30

# Размеры сетки
X_DIMENSION_GRID = 1.0E-3
Y_DIMENSION_GRID = 1.0E-2
Z_DIMENSION_GRID = 1.0E-2

# Число шагов по сетке для крупных частиц
X_STEP_NUMBER_GRID = 100
Y_STEP_NUMBER_GRID = 100
Z_STEP_NUMBER_GRID = 100

# Шаг по сетке
X_STEP = X_DIMENSION_GRID / X_STEP_NUMBER_GRID
Y_STEP = Y_DIMENSION_GRID / Y_STEP_NUMBER_GRID
Z_STEP = Z_DIMENSION_GRID / Z_STEP_NUMBER_GRID
STEPS = (X_STEP, Y_STEP, Z_STEP)

# Начальное значение времени метода установления
FAKE_TIME_STEP = 1.75E-9
while (sum(map(lambda step: FAKE_TIME_STEP/step**2, STEPS)) > 0.5):
    FAKE_TIME_STEP /= 2


# Начальное значение потенциала
POTENTIAL_BOUND_VALUE = 140

# Коэффициенты скорости для перехода к безразмерным величинам
ELECTRON_SPEED_DIMENSIONLESS_UNIT = 5.580192613E+5
CARBON_SPEED_DIMENSIONLESS_UNIT = 3760.714995
HELLIUM_SPEED_DIMENSIONLESS_UNIT = 1529.519939

# Коэффициенты для перехода к безразмерному пространству
SPACE_DIMENSIONLESS_UNIT = 2.211858553E-8

# Коэффициенты для перехода к безразмерной напряженности
INTENSITY_DIMENSIONLESS_UNIT = 4.002386595E+7

# Метод установления (epsilon)
ESTABLISHING_METHOD_ACCURACY = 0.01

# TODO Что это? Убрать?
Mt = 5.211008445E-12

# Сплоченность мелких частиц внутри крупной
PARTICLE_COHESION = 2.0  # 0..1

names = list(filter(lambda x: x.upper() == x, dir()))
CONSTANT_VALUES = dict(filter(lambda x: x[0] in names, globals().items()))
