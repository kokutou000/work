import numpy as np
import matplotlib.pyplot as plt

nnx = 514
nny = 514
head = ("head","<i4")
tail = ("tail","<i4")
sizebin2d = nnx*nny
chrsizebin2d = str(sizebin2d)
raw = ("data",chrsizebin2d+"<f8")
dt = np.dtype([head,raw,tail])
#
fdvx = open("vxca",'r')
chunkvx = np.fromfile(fdvx,dtype=dt,count=1)
datavx = chunkvx[0]["data"].reshape((nnx,nny),order="F")
fdvx.close()
#
fdvy = open("vyca",'r')
chunkvy = np.fromfile(fdvy,dtype=dt,count=1)
datavy = chunkvy[0]["data"].reshape((nnx,nny),order="F")
fdvy.close()
#
fdvx2d = open("vx2d","r")
chunkvx2d = np.fromfile(fdvx2d,dtype=dt,count=1)
datavx2d = chunkvx2d[0]["data"].reshape((nnx,nny),order="F")
fdvx2d.close()
#
plt.pcolormesh(datavx.T)
plt.colorbar()
plt.show()
#
plt.clf()
plt.plot(datavx[256,256:400])
plt.plot(datavx2d[256,256:400])
plt.show()
