import numpy as np
import matplotlib.pyplot as plt
import sys

nnx = 514
nny = 514
nnz = 257
head = ("head","<i4")
tail = ("tail","<i4")
sizebin3d = nnx*nny*nnz
sizebin2d = nnx*nny
chrsizebin3d = str(sizebin3d)
chrsizebin2d = str(sizebin2d)
raw3dbx = ("valbx",chrsizebin3d+"f8")
raw3dby = ("valby",chrsizebin3d+"f8")
raw3dbz = ("valbz",chrsizebin3d+"f8")
raw2dvx = ("valvx",chrsizebin2d+"f8")
raw2dvy = ("valvy",chrsizebin2d+"f8")
raw2dvz = ("valvz",chrsizebin2d+"f8")
dt3d = np.dtype([head,raw3dbx,raw3dby,raw3dbz,
                 raw2dvx,raw2dvy,raw2dvz,tail])
xc = np.genfromtxt("../coord.xgc")
yc = np.genfromtxt("../coord.ygc")

fddata = open("B3D_init","r")
chunkdata = np.fromfile(fddata,dtype=dt3d,count=1)
valbx = chunkdata[0]["valbx"].reshape((nnx,nny,nnz),order="F")
valby = chunkdata[0]["valby"].reshape((nnx,nny,nnz),order="F")
valbz = chunkdata[0]["valbz"].reshape((nnx,nny,nnz),order="F")
valvx = chunkdata[0]["valvx"].reshape((nnx,nny),order="F")
valvy = chunkdata[0]["valvy"].reshape((nnx,nny),order="F")
valvz = chunkdata[0]["valvz"].reshape((nnx,nny),order="F")
print(valbx.shape)
print(valvx.shape)

pltval3d = valbx
pltval2d = valvx
kz = 10
XC, YC = np.meshgrid(xc,yc)
# 3d plot
#plt.pcolormesh(pltval3d[:,:,kz].T,cmap="bwr")
plt.pcolormesh(XC,YC,pltval3d[:,:,kz].T,cmap="bwr")
plt.colorbar()
plt.show()
print(pltval3d[:5,300,kz])
print(pltval3d[508:,300,kz])
plt.close("all")
plt.plot(xc,pltval3d[:,300,kz],marker="+")
plt.plot(yc,pltval3d[300,:,kz],marker="^")
plt.show()
# 2d plot
# plt.pcolormesh(XC,YC,pltval2d[:,:].T,cmap="bwr")
# plt.colorbar()
# plt.xlim([-8.0,8.0])
# plt.ylim([-8.0,8.0])
# plt.show()

####################################################
sys.exit()
####################################################
tmpdbzdx = np.genfromtxt("dbzdx.dat")
tmpdbzdy = np.genfromtxt("dbzdy.dat")
tmpbz0 = np.genfromtxt("../../Bz_init.dat")
tmpdbzdytanh = np.genfromtxt("dbzdytanhbz.dat")
tmpdbzdycal = np.genfromtxt("dbzdycal.dat")
tmpvxtmp = np.genfromtxt("vxtmp.dat")
tmpvytmp = np.genfromtxt("vytmp.dat")

print(tmpdbzdx.shape)
print(nnx*nny)
dbzdx = tmpdbzdx[:,2].reshape(nnx,nny)
dbzdy = tmpdbzdy[:,2].reshape(nnx,nny)
bz0 = tmpbz0[:,2].reshape(nnx,nny)
dbzdytanh = tmpdbzdytanh[:,2].reshape(nnx,nny)
dbzdycal = tmpdbzdycal[:,1].reshape(nnx,nny)
vxtmp = tmpvxtmp.reshape(nnx,nny)
vytmp = tmpvytmp.reshape(nnx,nny)
print(dbzdx.shape)

xc = np.genfromtxt("../../coord.xgc")
yc = np.genfromtxt("../../coord.ygc")
XC, YC = np.meshgrid(xc,yc)

vvnorm = np.max(np.sqrt(vxtmp**2+vytmp**2))
print(vvnorm)

plt.close("all")
plt.pcolormesh(XC,YC,vxtmp.T,cmap="jet")
plt.colorbar()
plt.xlim([-5.0,5.0])
plt.ylim([-5.0,5.0])
plt.savefig("vxmap_modified.png")

sys.exit()

plt.close("all")
plt.pcolormesh(XC,YC,(vytmp/vvnorm).T,cmap="jet")
plt.colorbar()
plt.xlim([-5.0,5.0])
plt.ylim([-5.0,5.0])
plt.savefig("vymap_modified.png")

plt.close("all")
plt.plot(yc,bz0[257,:],marker="+")
plt.plot(yc,dbzdy[257,:]/(np.max(dbzdy[257,:])),marker="+")
plt.text(0.0,1.1,"max(dbzdy):"+str(np.max(dbzdy[257,:])))
plt.xlim([0.0,5.0])
plt.ylim([-0.5,1.0])
plt.savefig("bz_dbzdy__vs__y.png")

plt.close("all")
plt.plot(yc,bz0[257,:],marker="+")
plt.plot(yc,dbzdy[257,:]/(np.max(dbzdy[257,:])),marker="+")
plt.text(0.0,1.1,"max(dbzdy):"+str(np.max(dbzdy[257,:])))
plt.plot(yc,dbzdytanh[257,:]/(np.max(dbzdytanh[257,:])),marker="^")
plt.text(0.0,1.05,"max(dbzdytanh):"+str(np.max(dbzdytanh[257,:])))
plt.xlim([0.0,5.0])
plt.ylim([-0.5,1.0])
plt.savefig("dbzdytanh__vs__y.png")


plt.close("all")
plt.plot(yc,bz0[257,:],marker="+")
plt.plot(yc,dbzdy[257,:]/(np.max(dbzdy[257,:])),marker="+")
plt.text(0.0,1.1,"max(dbzdy):"+str(np.max(dbzdy[257,:])))
plt.plot(yc,dbzdytanh[257,:]/(np.max(dbzdytanh[257,:])),marker="^")
plt.text(0.0,1.05,"max(dbzdytanh):"+str(np.max(dbzdytanh[257,:])))
plt.plot(yc,vxtmp[257,:]/vvnorm,marker="+")
plt.xlim([0.0,5.0])
plt.ylim([-0.5,1.0])
plt.savefig("vxs__vs__y.png")
