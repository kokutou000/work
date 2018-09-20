import numpy as np
import matplotlib.pyplot as plt
import sys

nnx = 514
nny = 514
nnz = 257
head = ("head","<i4")
tail = ("tail","<i4")
sizebin3d = nnx*nny*nnz
chrsizebin3d = str(sizebin3d)
raw3d = ("val",chrsizebin3d+"f8")
dt3d = np.dtype([head,raw3d,tail])
xc = np.genfromtxt("../coord.xgc")
yc = np.genfromtxt("../coord.ygc")

fdbx = open("../Bz8bin_3d_init.dat","r")
chunkbx = np.fromfile(fdbx,dtype=dt3d,count=1)
valbx = chunkbx[0]["val"].reshape((nnx,nny,nnz),order="F")
print(valbx.shape)
kz = 0

plt.plot(yc,valbx[257,:,0],marker="+")
plt.show()
print(valbx[257,257:300,0])

sys.exit()
plt.pcolormesh(valbx[:,:,kz].T,cmap="bwr")
plt.colorbar()
plt.show()
print(valbx[:5,300,kz])
print(valbx[508:,300,kz])
plt.close("all")
plt.plot(xc,valbx[:,300,kz],marker="+")
plt.plot(yc,valbx[300,:,kz],marker="^")
plt.show()
plt.close("all")
plt.plot(xc,valbx[:,0,kz],marker="+",label="x")
plt.plot(yc,valbx[0,:,kz],marker="+",label="y")
plt.legend()
plt.show()
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
