import numpy as np
import matplotlib.pyplot as plt

rawxp = np.genfromtxt("../../../../3dbin/coord.xgc")
rawyp = np.genfromtxt("../../../../3dbin/coord.ygc")
vma=1.5; vmi=-1.5
infos = np.loadtxt("../../../../code_output/info",
                   dtype="string",comments="!")
tsini = int(infos[0]); tsfin = int(infos[1])
nnx = 1 + int(infos[3]) - int(infos[2])
nny = 1 + int(infos[5]) - int(infos[4])
nnz = 1 + int(infos[7]) - int(infos[6])
zcut = 0
###
var = "Tw2dbin."    
head = ("head","<i4"); tail = ("tail","<i4")
sizebin = nnx*nny; chrsizebin = str(sizebin)
raw = ("data",chrsizebin+"<f4"); dt = np.dtype([head,raw,tail])
###

for i in range(tsini,tsfin):
    ts = i
    cts = str(ts).zfill(3)

    fd = open("../../"+var+cts,'r')
    chunk = np.fromfile(fd,dtype=dt,count=1)
    data = chunk[0]["data"].reshape((nnx,nny),order="F")
    
    XP,YP = np.meshgrid(rawxp,rawyp)
    print("maxtw:", np.max(data[:,:]), "mintw:", np.min(data[:,:]))
    
    plt.clf()
    plt.pcolormesh(XP,YP,data[:,:].T,vmax=vma,vmin=vmi,cmap="coolwarm")
    plt.colorbar()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Tw dist. on xy plane @ z ="+str(zcut))
    plt.xlim([-4.0,4.0])
    plt.ylim([-4.0,4.0])
    plt.axes().set_aspect("equal","datalim")
    plt.savefig("tw2d"+cts,dpi=500)
