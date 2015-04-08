from constant import *
import numpy as np
import matplotlib.pyplot as plt

from constant import *

x = 0
normal_1 = np.random.normal(loc=x+X_STEP/2, scale=(X_STEP/6)**2, size=100)
normal_2 = np.random.normal(loc=x+X_STEP/2, scale=(X_STEP/6)**2, size=100)

plt.hist(normal_1, alpha=0.7, label='begin')
plt.hist(normal_2, alpha=0.7, label='end')
plt.legend()

plt.show()