import numpy as np
import matplotlib.pyplot as plt
import sys
import glob as gb
import os

rawxp = np.genfromtxt("../../../../3dbin/coord.xgc")
rawyp = np.genfromtxt("../../../../3dbin/coord.ygc")
rawzp = np.genfromtxt("../../../../3dbin/coord.zgc")
infos = np.loadtxt("../../../../code_output/info",
                   dtype="string",comments="!")
tsini = int(infos[0]); tsfin = int(infos[1])
nnx = 1 + int(infos[3]) - int(infos[2])
nny = 1 + int(infos[5]) - int(infos[4])
nnz = 1 + int(infos[7]) - int(infos[6])
nny0 = 513
nnz0 = 257
n_yz = np.genfromtxt("../../../../../../decay_index/n_yz.dat")
n_yz = n_yz.reshape(nny0+1,nnz0)
n_yz = n_yz[int(infos[2]):1+int(infos[3]),:int(infos[7])+1]
vma=20.0
vmi=0.0

sys.path.append(os.path.join(os.path.dirname(__file__),
                             "../../../../code_output"))
strtime = gb.glob("../../../../../TIME*")
rawtime = np.genfromtxt(strtime[0])
import get_info_restart as gir
timeinit=gir.get_info_restartfile("../../../"+infos[8],"LITTLE_ENDIAN")
timefin=gir.get_info_restartfile("../../../../../","BIG_ENDIAN")
timearr = rawtime[:,2]
timearr = np.append(timeinit,timearr)
timearr = np.append(timearr,timefin)
###
var = "JJ3dbin."
head = ("head","<i4");tail = ("tail","<i4")
sizebin = nnx*nny*nnz;chrsizebin = str(sizebin)
raw = ("data",chrsizebin+"<f4");dt = np.dtype([head,raw,tail])
###

for i in range(tsini,tsfin+1):
    ts = i
    cts = str(ts).zfill(3)
    print(cts)
        
    fd = open("../../"+var+cts,'r')
    chunk = np.fromfile(fd,dtype=dt,count=1)
    data = chunk[0]["data"].reshape((nnx,nny,nnz),order="F")
    
    YP, ZP = np.meshgrid(rawyp,rawzp)
#    sys.exit()

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    im = ax.pcolormesh(YP,ZP,data[(nnx-1)/2,:,:].T,vmax=vma,vmin=vmi,
                   cmap="jet")
    plt.colorbar(im, orientation="horizontal")
## draw decay index lines
    im2 = ax.contour(YP,ZP,n_yz.T,
                     levels=[0.5,1.0,1.5],
                     linestyles="dashed",
                     colors=["orange","yellow","red"])
    ax.set_xlim([-2.0,2.0])
    ax.set_ylim([0.0,2.0])
    ax.set_xlabel("y")
    ax.set_ylabel("z")
    ax.set_aspect("equal")
    ax.text(-2.0,-0.30,"t = "+str(timearr[ts-tsini])+" t_A")
    ax.text(-2.0,2.3, "max|J| = "+str(np.max(data)))
    plt.savefig("jjyz2d"+cts,dpi=500)

    plt.close()
