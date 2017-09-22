import numpy as np
import matplotlib.pyplot as plt
import matplotlib

a = np.loadtxt('sens_test.txt')
plt.scatter(a[:,0], a[:,1], 20, a[:,3])
plt.colorbar()
plt.show()


