# -*- coding: utf-8 -*-

import scipy
import scipy.stats
import numpy as np
import time


class MWranomizer(scipy.stats.rv_continuous):

    def __init__(self):
        ''' Maxwell distribution randomizer '''
        scipy.stats.rv_continuous.__init__(self)
        self.a = 0

    def _pdf(self, v, n, m, T):
        Boltzmann_constant = 1.380648813E-23
        k = Boltzmann_constant
        return n*(m/(2*np.pi*k*T))**(3/2)*np.exp(-m*v*v/(2*k*T))*4*np.pi*v**2

massCarbon = 12.011 * 1.67E-27 - 9.11E-31
temperature = 4200

maxwell_randomizer = MWranomizer()


size = 100000

start_time = time.time()

for i in range(size):
    speed = maxwell_randomizer.rvs(n=1, m=massCarbon, T=temperature, size=1)[0]
    print(speed)

elapsed_time = time.time() - start_time

start_time2 = time.time()
speed = maxwell_randomizer.rvs(n=1, m=massCarbon, T=temperature, size=size)
for i in range(size):
    print(speed[i])
elapsed_time2 = time.time() - start_time2



print(elapsed_time, elapsed_time2)
