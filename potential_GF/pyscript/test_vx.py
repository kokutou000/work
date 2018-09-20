import numpy as np
import matplotlib.pyplot as plt

rawxp = np.genfromtxt("../coord.xgc")
rawyp = np.genfromtxt("../coord.ygc")
rawzp = np.genfromtxt("../coord.zgc")
rawvx0 = np.genfromtxt("../Vx_init.dat")

nnx = 514
nny = 514
vx_init = rawvx0[:,2].reshape(nnx,nny)

print(vx_init.shape)

XP_z, YP_z = np.meshgrid(rawxp,rawyp)

plt.pcolormesh(XP_z, YP_z, vx_init.T,
               cmap="jet")
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.title("Vx dist. on z=0")
plt.xlim([-5.0,5.0])
plt.ylim([-5.0,5.0])
#plt.show()
plt.savefig("vxmap.png",dpi=500)
