import numpy as np
import matplotlib.pyplot as plt
import sys
import glob as gb
import os

#rawxp = np.genfromtxt("../../../../binfile/coord.xgc")
#rawyp = np.genfromtxt("../../../../binfile/coord.ygc")
rawxp = np.genfromtxt("../../../../3dbin/coord.xgc")
rawyp = np.genfromtxt("../../../../3dbin/coord.ygc")
rawzp = np.genfromtxt("../../../../3dbin/coord.zgc")
vma=1.5
vmi=-1.5
infos = np.genfromtxt("../../../../code_output/info",
                      dtype="string",
                      comments="!")
sys.path.append(os.path.join(os.path.dirname(__file__),
                             "../../../../code_output"))
tsini=int(infos[0]); tsfin=int(infos[1])
height_axis = np.zeros(tsfin-tsini+1)
strtime=gb.glob("../../../../../TIME*")
rawtime=np.genfromtxt(strtime[0])
time=rawtime[:,2]
import get_info_restart as gir
timeinit=gir.get_info_restartfile("../../../"+infos[8],"LITTLE_ENDIAN")
timefin=gir.get_info_restartfile("../../../../../","BIG_ENDIAN")
time = np.append(timeinit,time)
time = np.append(time,timefin)
##### get information of z(n=1.0) and z(n=1.5) #####
n_at0 = np.genfromtxt("../../../../../../decay_index/n_at_center.dat")
zn10 = np.argmin(np.abs(n_at0[:,1]-1.0))
zn15 = np.argmin(np.abs(n_at0[:,1]-1.5))
####################################################
nnx = int(infos[3]) - int(infos[2]) + 1
nny = int(infos[5]) - int(infos[4]) + 1
nnz = int(infos[7]) - int(infos[6]) + 1
nnx2 = int(nnx/2); nny2 = int(nny/2)
zcut = 0
var = "Bybin_3d_R."
head = ("head","<i4"); tail = ("tail","<i4")
sizebin = nnx*nny*nnz; chrsizebin = str(sizebin)
raw = ("data",chrsizebin+"<f4"); dt = np.dtype([head,raw,tail])

for i in range(tsini,tsfin+1):
    ts = i
    cts = str(ts).zfill(3)
    
    fd = open("../../../../3dbin/"+var+cts,'r')
    chunk = np.fromfile(fd,dtype=dt,count=1)
    data = chunk[0]["data"].reshape((nnx,nny,nnz),order="F")
    
    by_z = data[nnx2,nny2,:]
    tmpkz = 0
    for kz in range(1,nnz-1):
        if (by_z[kz+1]*by_z[kz-1])<0.0:
            tmpkz = kz
#    tmp_height = rawzp[np.argmin(np.abs(data[nnx2,nny2,:]))]
    tmp_height = rawzp[tmpkz]
    print(ts,tmp_height)
    height_axis[i-tsini] = tmp_height
    
endt_inj = float(infos[10])
plt.clf()
### linear plot
plt.plot(time,height_axis,marker="+")
plt.axhline(n_at0[zn10,0],linestyle="dashed",color="yellow")
plt.axhline(n_at0[zn15,0],linestyle="dashed",color="red")
plt.axvline(endt_inj,linestyle="dashed",color="black")
plt.xlabel("time")
plt.ylabel("height of FR axis")
plt.ylim([0.0,2.0])
plt.savefig("evol_FR_axis.png",dpi=500)
###
plt.close("all")
### log10 plot
plt.plot(time,np.log10(height_axis),marker="+")
plt.axvline(endt_inj,linestyle="dashed",color="black")
plt.xlabel("time")
plt.ylabel("height of FR axis")
plt.savefig("log10_evol_FR_axis.png",dpi=500)

np.savetxt("data_height.txt",height_axis)
