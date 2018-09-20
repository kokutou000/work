import numpy as np
import matplotlib.pyplot as plt

rawxp = np.genfromtxt("../coord.xgc")
rawyp = np.genfromtxt("../coord.ygc")
rawzp = np.genfromtxt("../coord.zgc")
rawvx0 = np.genfromtxt("../Vx_init.dat")
rawbz0 = np.genfromtxt("../Bz_init.dat")

nnx = 514
nny = 514
vx_init = rawvx0[:,2].reshape(nnx,nny)
bz_init = rawbz0[:,2].reshape(nnx,nny)

print(vx_init.shape)

XP_z, YP_z = np.meshgrid(rawxp,rawyp)

plt.plot(rawyp, vx_init[nnx/2,:]*100.0,marker="^",label="Vx*100")
plt.plot(rawyp, bz_init[nnx/2,:],marker=".",label="Bz")
plt.xlim([0.0,4.0])
plt.xlabel("y")
plt.ylabel("Bz or Vx*100")
plt.legend(loc="lower right")
#plt.show()
plt.savefig("plot1d_bzvx.png")

