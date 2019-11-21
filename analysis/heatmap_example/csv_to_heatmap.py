import numpy as np
#ifrom mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
#from scipy.interpolate import griddata
from matplotlib import cm
import matplotlib.tri as tri

# Load data from CSV
dat = np.genfromtxt('data_scan.dat', delimiter=', ',skip_header=0)
X_dat = dat[:,0]
Y_dat = dat[:,1]
Z_dat = dat[:,2]

# Convert from pandas dataframes to numpy arrays
X, Y, Z, = np.array([]), np.array([]), np.array([])
for i in range(len(X_dat)):
        X = np.append(X,X_dat[i])
        Y = np.append(Y,Y_dat[i])
        Z = np.append(Z,Z_dat[i])

triang = tri.Triangulation(X, Y)

#plt.gca().set_aspect('equal')
plt.tricontourf(triang, Z, 1000, cmap=cm.plasma)
plt.colorbar()
#plt.tricontour(triang, Z, colors='k')
plt.xlim(6.1, 6.25)
plt.ylim(90, 95)

plt.show()

# create x-y points to be used in heatmap
#xi = np.linspace(X.min(),X.max(),1000)
#yi = np.linspace(Y.min(),Y.max(),1000)

# Z is a matrix of x-y values
#zi = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method='cubic')

#print xi
#print yi
#print zi

#print np.isnan(zi).sum()


# I control the range of my colorbar by removing data 
# outside of my range of interest
#zmin = -105
#zmax = -100
#zi[(zi<zmin) | (zi>zmax)] = None

# Create the contour plot
#CS = plt.contourf(xi, yi, zi, 15, cmap=plt.cm.rainbow,
#                  vmax=zmax, vmin=zmin)

#fig = plt.figure()
#ax = fig.gca(projection='3d')
# Plot the surface.
#surf = ax.plot_surface(X, Y, Z)
#plt.show()

#xi, yi = np.meshgrid(np.linspace(X.min(),X.max(),1000),np.linspace(Y.min(),Y.max(),1000))
#interp_cubic_geom = tri.CubicTriInterpolator(triang, Z, kind='geom')
#zi_cubic_geom = interp_cubic_geom(xi, yi)

#plt.contour(xi, yi, zi_cubic_geom)
#plt.plot(xi, yi, 'k-', lw=0.5, alpha=0.5)
#plt.plot(xi.T, yi.T, 'k-', lw=0.5, alpha=0.5)


#fig = plt.figure()
#ax = fig.gca(projection='3d')
#surf = ax.plot_surface(xi, yi, zi_cubic_geom, linewidth=0)
