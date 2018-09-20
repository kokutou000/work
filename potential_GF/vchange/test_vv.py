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
#dirname = "nearPIL1"
dirname = "nomask"

fdvx = open("../"+dirname+"/vxca","r")
chunkvx = np.fromfile(fdvx,dtype=dt2d,count=1)
valvx = chunkvx[0]["val"].reshape((nx,ny),order="F")
#
fdvx0 = open("../"+dirname+"/vx2d","r")
chunkvx0 = np.fromfile(fdvx0,dtype=dt2d,count=1)
valvx0 = chunkvx0[0]["val"].reshape((nx,ny),order="F")
#
fdvy = open("../"+dirname+"/vyca","r")
chunkvy = np.fromfile(fdvy,dtype=dt2d,count=1)
valvy = chunkvy[0]["val"].reshape((nx,ny),order="F")
#
fdvy0 = open("../"+dirname+"/vy2d","r")
chunkvy0 = np.fromfile(fdvy0,dtype=dt2d,count=1)
valvy0 = chunkvy0[0]["val"].reshape((nx,ny),order="F")

valvv0 = (valvx0**2 + valvy0**2)**0.5
valvv = (valvx**2 + valvy**2)**0.5

print(valvx.shape)
rawxp = np.genfromtxt("../../input/coord.xgc")
rawyp = np.genfromtxt("../../input/coord.ygc")
print(rawxp)

#plt.clf()
fig = plt.figure(figsize=(12,6))
axarr0 = fig.add_subplot(121)
im0 = plt.pcolormesh(rawxp,rawyp,valvx0.T,cmap="jet")
plt.colorbar()
# plt.contour(rawxp,rawyp,
#             valvv0.T,
#             levels=[0.0,0.0025,0.005,0.0075],
#             linestyles="dashed")
axarr1 = fig.add_subplot(122)
im1 = plt.pcolormesh(rawxp,rawyp,valvx.T,cmap="jet")
plt.colorbar()
# plt.contour(rawxp,rawyp,
#             valvv.T,
#             levels=[0.0,0.0025,0.005,0.0075],
#             linestyles="dashed")
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
fig = plt.figure(figsize=(12,6))
axarr0 = fig.add_subplot(121)
im0 = plt.pcolormesh(rawxp,rawyp,valvy0.T,cmap="jet")
plt.colorbar()
# plt.contour(rawxp,rawyp,
#             valvv0.T,
#             levels=[0.0,0.0025,0.005,0.0075],
#             linestyles="dashed")
axarr1 = fig.add_subplot(122)
im1 = plt.pcolormesh(rawxp,rawyp,valvy.T,cmap="jet")
plt.colorbar()
# plt.contour(rawxp,rawyp,
#             valvv.T,
#             levels=[0.0,0.0025,0.005,0.0075],
#             linestyles="dashed")
axarr0.set_xlim([-4,4])
axarr0.set_ylim([-4,4])
axarr0.set_aspect("equal")
axarr1.set_xlim([-4,4])
axarr1.set_ylim([-4,4])
axarr1.set_aspect("equal")
plt.savefig("../"+dirname+"/cmapvy.png",dpi=500)
#plt.show()
plt.close("all")
#
xcut = (nx-1)/2
plt.plot(rawyp,valvx[xcut,:],marker="^")
plt.plot(rawyp,valvx0[xcut,:],marker="+")
plt.xlim([0,4.0])
plt.savefig("../"+dirname+"/plotvx1d.png")
#plt.show()
