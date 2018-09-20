import numpy as np
import matplotlib.pyplot as plt

#nnx = 514
#nny = 514
#nnz = 257
nnx = 514
nny = 514
nnz = 257
head = ("head","<i4")
tail = ("tail","<i4")
sizebin3d = nnx*nny*nnz
chrsizebin3d = str(sizebin3d)
raw3dcx = ("data3dcx",chrsizebin3d+"f8")
raw3dcy = ("data3dcy",chrsizebin3d+"f8")
raw3dcz = ("data3dcz",chrsizebin3d+"f8")
dtcx = np.dtype([head,raw3dcx,tail])
dtcy = np.dtype([head,raw3dcy,tail])
dtcz = np.dtype([head,raw3dcz,tail])

fdcx = open("../cxd0_direct","r")
chunkcx = np.fromfile(fdcx,dtype=dtcx,count=1)
datacx = chunkcx[0]["data3dcx"].reshape((nnx,nny,nnz),order="F")
fdcy = open("../cyd0_direct","r")
chunkcy = np.fromfile(fdcy,dtype=dtcy,count=1)
datacy = chunkcy[0]["data3dcy"].reshape((nnx,nny,nnz),order="F")
fdcz = open("../czd0_direct","r")
chunkcz = np.fromfile(fdcz,dtype=dtcz,count=1)
datacz = chunkcz[0]["data3dcz"].reshape((nnx,nny,nnz),order="F")

rawxp = np.genfromtxt("../coord.xgc")
rawyp = np.genfromtxt("../coord.ygc")
rawzp = np.genfromtxt("../coord.zgc")
XP_z, YP_z = np.meshgrid(rawxp,rawyp)
YP_x, ZP_x = np.meshgrid(rawyp,rawzp)
XP_y, ZP_y = np.meshgrid(rawxp,rawzp)

print(datacx.shape)

#plt.pcolormesh(XP_z,YP_z,databz[:,:,2].T)
#plt.pcolormesh(XP_y,ZP_y,databz[:,64,:].T)
#plt.pcolormesh(YP_x,ZP_x,databz[64,:,:].T)
#plt.pcolormesh(YP_x,ZP_x,databy[64,:,:].T)
plt.pcolormesh(XP_z,YP_z,datacx[:,:,0].T)
plt.colorbar()
#plt.ylim([0,5])
plt.show()

print("max(cx):",np.max(datacx[:,:,0:5]), "min(cx):", np.min(datacx[:,:,0:5]))
print("max(cy):",np.max(datacy[:,:,0:5]), "min(cy):", np.min(datacy[:,:,0:5]))
print("max(cz):",np.max(datacz[:,:,0:5]), "min(cz):", np.min(datacz[:,:,0:5]))
