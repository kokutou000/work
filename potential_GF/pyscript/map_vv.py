import numpy as np
import matplotlib.pyplot as plt

rawxp = np.genfromtxt("../coord.xgc")
rawyp = np.genfromtxt("../coord.ygc")
rawzp = np.genfromtxt("../coord.zgc")
rawvx0 = np.genfromtxt("../Vx_init.dat")
rawvy0 = np.genfromtxt("../Vy_init.dat")
rawbz0 = np.genfromtxt("../Bz_init.dat")

nnx = 514
nny = 514
vx_init = rawvx0[:,2].reshape(nnx,nny)
vy_init = rawvy0[:,2].reshape(nnx,nny)
bz_init = rawbz0[:,2].reshape(nnx,nny)
vv_init = np.sqrt(vx_init**2 + vy_init**2)

print(vx_init.shape)

XP_z, YP_z = np.meshgrid(rawxp,rawyp)

plt.pcolormesh(XP_z, YP_z, vv_init.T,
               cmap="jet")
plt.colorbar()
plt.contour(XP_z, YP_z, bz_init.T,
            levels=[-0.8,-0.5,0.0,0.5,0.8],
            linestyles="dashed",
            linewidths=1.0)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Vx dist. on z=0")
plt.xlim([-5.0,5.0])
plt.ylim([-5.0,5.0])
plt.show()
