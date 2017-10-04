import numpy as np
import matplotlib.pyplot as plt
import matplotlib

a = np.loadtxt('gpts_sens_test.txt')
b = np.loadtxt('bpts_sens_test.txt')

plt.figure(1)
plt.scatter(a[:,0], a[:,1], 20, a[:,3])
plt.colorbar()

plt.figure(2)
plt.scatter(b[:,0], b[:,1], 20, b[:,2])
plt.colorbar()
plt.show()


