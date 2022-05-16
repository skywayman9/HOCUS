# from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
# from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate as interpolate


# # For latex font, i guess so 
plt.rc('text', usetex=True)
plt.rc('font', family='arial')

plt.rcParams.update({'font.size': 12})


# x,y,d,p,u,v=np.loadtxt('soln.txt', delimiter=None, unpack=True)

x,y,d,u,v,p=np.loadtxt('Rslt0008.tec', delimiter=None, unpack=True,skiprows=3)


g=200
k=400


x = x.reshape(g,k)
y = y.reshape(g,k)
d = d.reshape(g,k)
p = p.reshape(g,k)
u = u.reshape(g,k)
v = v.reshape(g,k)

# plt.xlim(0.0,2.0)
# plt.ylim(0.0,1.0)
plt.contour(x,y,d,24,linewidths=0.5,colors=('k'))

plt.ylabel(r'\textbf{y}')
plt.xlabel(r'\textbf{x}')
fig1 = plt.gcf()
# plt.colorbar()
fig1.set_size_inches(w=6,h=3, forward=True)
fig1.savefig('Shock_vortex.pdf', dpi=600,bbox_inches='tight', pad_inches = 0)
plt.show()
