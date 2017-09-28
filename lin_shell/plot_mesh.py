import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib

node = np.loadtxt('node.txt')
elem = np.loadtxt('elem.txt')

fig = plt.figure()
ax = fig.gca(projection='3d')
for ii in range(0, elem.shape[0]):
    elem_id = elem[ii,:]
    elem_nodes = np.array([int(elem_id[0]),int(elem_id[1]),int(elem_id[2]),int(elem_id[3]),int(elem_id[0])])
    ax.plot(node[elem_nodes,0],node[elem_nodes,1],node[elem_nodes,2],'k')

plt.show()

