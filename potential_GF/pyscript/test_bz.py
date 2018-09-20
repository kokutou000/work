import numpy as np
import matplotlib.pyplot as plt

rawxp = np.genfromtxt("../coord.xgc")
rawyp = np.genfromtxt("../coord.ygc")
rawzp = np.genfromtxt("../coord.zgc")
rawbz0 = np.genfromtxt("../Bz_init.dat")

nnx = 514
nny = 514
bz_init = rawbz0[:,2].reshape(nnx,nny)

print(bz_init.shape)

XP_z, YP_z = np.meshgrid(rawxp,rawyp)

plt.pcolormesh(XP_z, YP_z, bz_init.T,
               cmap="gray")
plt.colorbar()
plt.contour(XP_z,YP_z,bz_init.T,
            levels=[0.0,0.25,0.5,0.75],
            linestyles="dashed")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Bz dist. on z=0")
plt.xlim([-5.0,5.0])
plt.ylim([-5.0,5.0])
#plt.show(dpi=300)
plt.savefig("bz_init.png",dpi=500)

print(bz_init[nnx/2,(nny/2)-2:(nny/2)+3])
print(bz_init[nnx/2,nny/2+1],bz_init[nnx/2,nny/2-1])
print(rawyp[nny/2+1],rawyp[nny/2-1])
dbzdy = (bz_init[nnx/2,nny/2+1]-bz_init[nnx/2,nny/2-1])/(rawyp[nny/2+1]-rawyp[nny/2-1])
print(dbzdy)
