import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("newfunction.dat")
print(data.shape)
data2 = np.genfromtxt("dbzdycal.dat")

nx = 514
ny = 514
print(nx*ny)
newf = data[:,1].reshape(nx,ny)
dbzdycal = data2[:,1].reshape(nx,ny)

#plt.pcolormesh(newf.T,vmax=1.0)
plt.pcolormesh(dbzdycal.T,vmax=1.0)
plt.colorbar()
plt.show()
print(dbzdycal[257,279])
