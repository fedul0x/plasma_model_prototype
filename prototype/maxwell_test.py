# -*- coding: utf-8 -*-

import scipy
import time
import scipy.stats
from scipy.stats import maxwell
import numpy as np
import matplotlib.pyplot as plt

from constant import *


class MWranomizer(scipy.stats.rv_continuous):

    def __init__(self):
        ''' Maxwell distribution randomizer '''
        scipy.stats.rv_continuous.__init__(self)
        self.a = 0

    def _pdf(self, v):
        Boltzmann_constant = 1.380648813E-23
        k = Boltzmann_constant
        n = 1
        m = CARBONS_MASS
        T = TEMPERATURE
        return n*(m/(2*np.pi*k*T))**(3/2)*np.exp(-m*v*v/(2*k*T))*4*np.pi*v**2

massCarbon = 12.011 * 1.67E-27 - 9.11E-31
temperature = 4200

maxwell_randomizer = MWranomizer()


# size = 100000
size = 10

start_time = time.time()

for i in range(size):
    speed = maxwell_randomizer.rvs(size=1)[0]
    print(speed)

elapsed_time = time.time() - start_time

start_time2 = time.time()
speed_1 = maxwell_randomizer.rvs(size=size)/CARBON_SPEED_DIMENSIONLESS_UNIT
for i in range(size):
    print(speed_1[i])
elapsed_time2 = time.time() - start_time2

start_time3 = time.time()
for i in range(size):
    speed = maxwell.rvs(scale=CARBONS_MAX_SPEED/1.5957691216057308)
    print(speed)
elapsed_time3 = time.time() - start_time3

start_time4 = time.time()
# speed_2 = maxwell.rvs(scale=CARBONS_MAX_SPEED/1.5957691216057308, size=size)
speed_2 = maxwell.rvs(size=size)
for i in range(size):
    print(speed_2[i])
elapsed_time4 = time.time() - start_time4

print('======================')
for i in range(10):
    print('{} = {}'.format(maxwell_randomizer.stats(moments='mv'),
        maxwell.stats(scale=CARBONS_MAX_SPEED/1.5957691216057308, moments='mv')))

fig, (ax1, ax2) = plt.subplots(2)
ax1.hist(speed_1, normed=True, histtype='stepfilled')
ax2.hist(speed_2, normed=True, histtype='stepfilled')
plt.show()



print(elapsed_time, elapsed_time2, elapsed_time3, elapsed_time4)
