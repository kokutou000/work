import numpy as np
import matplotlib.pyplot as plt

rawxc = np.genfromtxt("../coord.xgc")
rawyc = np.genfromtxt("../coord.ygc")
rawzc = np.genfromtxt("../coord.zgc")

nnx = 514
nny = 514
nnz = 257
para_vmax=2.0
para_vmin=0.0

rawn_yz = np.genfromtxt("n_yz.dat")
nc_yz = rawn_yz.reshape(nny,nnz)
YP_x, ZP_x = np.meshgrid(rawyc,rawzc)

plt.pcolormesh(YP_x, ZP_x,nc_yz.T,
               cmap="Blues",
               vmax=para_vmax,
               vmin=para_vmin)
plt.colorbar()
plt.contour(YP_x, ZP_x, nc_yz.T,
            levels=[0.0,0.5,1.0,1.5],
            linestyles="--")
plt.xlim([-2.0,2.0])
plt.ylim([0.0,4.0])
plt.xlabel("y")
plt.ylabel("z")
plt.title("decay index on yz plane")
plt.savefig("decay_index_yz.png")
