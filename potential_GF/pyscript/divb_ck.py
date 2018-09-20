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
raw3ddbxdx = ("data3ddbxdx",chrsizebin3d+"f8")
raw3ddbydy = ("data3ddbydy",chrsizebin3d+"f8")
raw3ddbzdz = ("data3ddbzdz",chrsizebin3d+"f8")
raw3ddivb = ("data3ddivb",chrsizebin3d+"f8")
dtdbxdx = np.dtype([head,raw3ddbxdx,tail])
dtdbydy = np.dtype([head,raw3ddbydy,tail])
dtdbzdz = np.dtype([head,raw3ddbzdz,tail])
dtdivb = np.dtype([head,raw3ddivb,tail])

fddbxdx = open("../dbxdx_debug","r")
chunkdbxdx = np.fromfile(fddbxdx,dtype=dtdbxdx,count=1)
datadbxdx = chunkdbxdx[0]["data3ddbxdx"].reshape((nnx,nny,nnz),order="F")
fddbydy = open("../dbydy_debug","r")
chunkdbydy = np.fromfile(fddbydy,dtype=dtdbydy,count=1)
datadbydy = chunkdbydy[0]["data3ddbydy"].reshape((nnx,nny,nnz),order="F")
fddbzdz = open("../dbzdz_debug","r")
chunkdbzdz = np.fromfile(fddbzdz,dtype=dtdbzdz,count=1)
datadbzdz = chunkdbzdz[0]["data3ddbzdz"].reshape((nnx,nny,nnz),order="F")
fddivb = open("../divb_debug","r")
chunkdivb = np.fromfile(fddivb,dtype=dtdivb,count=1)
datadivb = chunkdivb[0]["data3ddivb"].reshape((nnx,nny,nnz),order="F")

rawxp = np.genfromtxt("../coord.xgc")
rawyp = np.genfromtxt("../coord.ygc")
rawzp = np.genfromtxt("../coord.zgc")
XP_z, YP_z = np.meshgrid(rawxp,rawyp)
YP_x, ZP_x = np.meshgrid(rawyp,rawzp)
XP_y, ZP_y = np.meshgrid(rawxp,rawzp)

#plt.pcolormesh(XP_z,YP_z,databz[:,:,2].T)
#plt.pcolormesh(XP_y,ZP_y,databz[:,64,:].T)
#plt.pcolormesh(YP_x,ZP_x,databz[64,:,:].T)
#plt.pcolormesh(YP_x,ZP_x,databy[64,:,:].T)
plt.pcolormesh(XP_z,YP_z,datadivb[:,:,3].T)
plt.colorbar()
#plt.ylim([0,5])
plt.show()

print("max(divb):",np.max(datadivb[:,:,0:5]),
      "min(divb):",np.min(datadivb[:,:,0:5]))
#print("max(cx):",np.max(datacx[:,:,0:5]), "min(cx):", np.min(datacx[:,:,0:5]))
#print("max(cy):",np.max(datacy[:,:,0:5]), "min(cy):", np.min(datacy[:,:,0:5]))
#print("max(cz):",np.max(datacz[:,:,0:5]), "min(cz):", np.min(datacz[:,:,0:5]))
