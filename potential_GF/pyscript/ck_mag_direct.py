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
rawbx4b = ("databx4b",chrsizebin3d+"f8")
rawby4b = ("databy4b",chrsizebin3d+"f8")
rawbz4b = ("databz4b",chrsizebin3d+"f8")
dtbx = np.dtype([head,rawbx4b,tail])
dtby = np.dtype([head,rawby4b,tail])
dtbz = np.dtype([head,rawbz4b,tail])

fdbx = open("../bxd0_direct","r")
chunkbx = np.fromfile(fdbx,dtype=dtbx,count=1)
databx = chunkbx[0]["databx4b"].reshape((nnx,nny,nnz),order="F")
fdby = open("../byd0_direct","r")
chunkby = np.fromfile(fdby,dtype=dtby,count=1)
databy = chunkby[0]["databy4b"].reshape((nnx,nny,nnz),order="F")
fdbz = open("../bzd0_direct","r")
chunkbz = np.fromfile(fdbz,dtype=dtbz,count=1)
databz = chunkbz[0]["databz4b"].reshape((nnx,nny,nnz),order="F")

rawxp = np.genfromtxt("../coord.xgc")
rawyp = np.genfromtxt("../coord.ygc")
rawzp = np.genfromtxt("../coord.zgc")
XP_z, YP_z = np.meshgrid(rawxp,rawyp)
YP_x, ZP_x = np.meshgrid(rawyp,rawzp)
XP_y, ZP_y = np.meshgrid(rawxp,rawzp)

databb = (databx**2 + databy**2 + databz**2)**(0.5)

#for jy in range(256,260):
#    plt.pcolormesh(XP_y,ZP_y,databz[:,jy,:].T)
for kz in range(0,3):
    plt.pcolormesh(XP_z, YP_z, databz[:,:,kz].T)
    plt.colorbar()
    plt.ylim([-4.0,4.0])
    plt.xlim([-4.0,4.0])
#    plt.ylim([-1.0,1.0])
    plt.show()
    plt.close("all")


print("max(bx):",np.max(databx[:,:,0:5]), "min(bx):", np.min(databx[:,:,0:5]))
print("max(by):",np.max(databy[:,:,0:5]), "min(by):", np.min(databy[:,:,0:5]))
print("max(bz):",np.max(databz[:,:,0:5]), "min(bz):", np.min(databz[:,:,0:5]))
