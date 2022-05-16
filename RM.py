# from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

# # For latex font, i guess so 
plt.rc('text', usetex=True)
plt.rc('font', family='arial')
#Set global matplotlib parameters in script or in /home/$USER/.matplotlib/matplotlibrc
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams.update({'font.size': 12})
plt.rcParams['figure.figsize'] = 12, 4

# x,y,d,p,u,v=np.loadtxt('soln.txt', delimiter=None, unpack=True)

x,y,d,u,v,p=np.loadtxt('Rslt0031.tec', delimiter=None, unpack=True,skiprows=3)


g=80
k=320


x = x.reshape(g,k)
y = y.reshape(g,k)
d = d.reshape(g,k)
p = p.reshape(g,k)
u = u.reshape(g,k)
v = v.reshape(g,k)

plt.xlim(0.0,1.5)
plt.ylim(0.0,1.0)

plt.contour(x,y,d,24,linewidths=0.5,colors=('k'))

plt.ylabel(r'\textbf{y}')
plt.xlabel(r'\textbf{x}')
fig1 = plt.gcf()
# plt.colorbar()
fig1.set_size_inches(w=3,h=3/1.5, forward=True)
fig1.savefig('Meshkov.pdf', dpi=600,bbox_inches='tight', pad_inches = 0)
plt.show()
