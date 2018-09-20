import numpy as np
import matplotlib.pyplot as plt

nx = 514
ny = 514
nz = 257
head = ("head","<i4")
tail = ("tail","<i4")
sizebin2d = nx*ny
chrsizebin2d = str(sizebin2d)
raw2d = ("val",chrsizebin2d+"f8")
dt2d = np.dtype([head,raw2d,tail])
dirname = "nearPIL1.1"
#dirname = "notanh"

fdvx = open("../"+dirname+"/vxca","r")
chunkvx = np.fromfile(fdvx,dtype=dt2d,count=1)
valvx = chunkvx[0]["val"].reshape((nx,ny),order="F")
#
fdvx0 = open("../"+dirname+"/vx2d","r")
chunkvx0 = np.fromfile(fdvx0,dtype=dt2d,count=1)
valvx0 = chunkvx0[0]["val"].reshape((nx,ny),order="F")

print(valvx.shape)
rawxp = np.genfromtxt("../../input/coord.xgc")
rawyp = np.genfromtxt("../../input/coord.ygc")
print(rawxp)

#plt.clf()
fig = plt.figure(figsize=(12,6))
axarr0 = fig.add_subplot(121)
im0 = plt.pcolormesh(rawxp,rawyp,valvx0.T,cmap="jet")
plt.colorbar()
axarr1 = fig.add_subplot(122)
im1 = plt.pcolormesh(rawxp,rawyp,valvx.T,cmap="jet")
plt.colorbar()
axarr0.set_xlim([-4,4])
axarr0.set_ylim([-4,4])
axarr0.set_aspect("equal")
axarr1.set_xlim([-4,4])
axarr1.set_ylim([-4,4])
axarr1.set_aspect("equal")
#plt.show()
plt.savefig("../"+dirname+"/cmapvx.png",dpi=500)
plt.close("all")
#
xcut = (nx-1)/2
plt.plot(rawyp,valvx[xcut,:],marker="^")
plt.plot(rawyp,valvx0[xcut,:],marker="+")
plt.xlim([0,4.0])
plt.savefig("../"+dirname+"/plotvx1d.png")
#plt.show()
