import numpy as np
import matplotlib.pyplot as plt

rawxc = np.genfromtxt("../coord.xgc")
rawzc = np.genfromtxt("../coord.zgc")

nnx = 514
nnz = 257
para_vmax=2.0
para_vmin=0.0

rawn = np.genfromtxt("n_zx.dat")
nc = rawn.reshape(nnx,nnz)
XP_y, ZP_y = np.meshgrid(rawxc,rawzc)

plt.pcolormesh(XP_y, ZP_y,nc.T,
               cmap="Blues",
               vmax=para_vmax,
               vmin=para_vmin)
plt.colorbar()
plt.contour(XP_y, ZP_y, nc.T,
            levels=[0.0,0.5,1.0,1.5],
            linestyles="--")
plt.xlim([-2.0,2.0])
plt.ylim([0.0,4.0])
plt.xlabel("x along PIL")
plt.ylabel("z")
plt.title("decay index")
plt.savefig("decay_index_zx.png")
