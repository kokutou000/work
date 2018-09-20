import numpy as np
import matplotlib.pyplot as plt

rawxp = np.genfromtxt("../coord.xgc")
rawyp = np.genfromtxt("../coord.ygc")
rawzp = np.genfromtxt("../coord.zgc")
rawvy0 = np.genfromtxt("../Vy_init.dat")

nnx = 514
nny = 514
vy_init = rawvy0[:,2].reshape(nnx,nny)

print(vy_init.shape)

XP_z, YP_z = np.meshgrid(rawxp,rawyp)

plt.pcolormesh(XP_z, YP_z, vy_init.T,
               cmap="jet")
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.title("Vy dist. on z=0")
plt.xlim([-5.0,5.0])
plt.ylim([-5.0,5.0])
#plt.show()
plt.savefig("vymap.png",dpi=500)
