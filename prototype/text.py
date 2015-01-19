import scipy.stats
import matplotlib.pyplot as plt
import numpy as np

class MWranomizer(scipy.stats.rv_continuous):
    def __init__(self):
        scipy.stats.rv_continuous.__init__(self)
        self.a = 0

    def _pdf(self, v, n, m, T):
        Boltzmann_constant = 1.380648813E-23
        k = Boltzmann_constant
        return n*(m/(2*np.pi*k*T))**(3/2)*np.exp(-m*v*v/(2*k*T))*4*np.pi*v**2
        # return n*np.exp(v)*4*np.pi*v**2


m = MWranomizer()
# m.a = 0
# h = m.rvs(n=1000, m=1E-17, T=273, size=1000)
h = m.rvs(n=1, m=2.005745900E-26, T=4200, size=5000)
# print(h)
plt.hist(h)
plt.show()

