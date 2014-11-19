import matplotlib.pylab as plt
import numpy as np

Ex = np.load('outs/Fields1_Ex_0.356242801958_x0.382344502427_y.npy')
Ey = np.load('outs/Fields1_Ey_0.356242801958_x0.382344502427_y.npy')

Irrad = np.abs(Ex)**2 + np.abs(Ey)**2

fig = plt.figure()
ax = fig.add_subplot(111, aspect = 'equal')
ax.pcolor(Irrad)
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("field_dist.png")
fig.show()
