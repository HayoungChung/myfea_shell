import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# for ii in range(1,5):
ii = int(19)
a = np.loadtxt('step_%d.gptSens.txt'%ii)
b = np.loadtxt('step_%d.bptsSens.txt'%ii)

plt.figure(1)
plt.scatter(a[:,0], a[:,1], 20, a[:,3])
plt.colorbar()

plt.figure(2)
plt.scatter(b[:,0], b[:,1], 20, b[:,2])
plt.colorbar()
plt.show()
