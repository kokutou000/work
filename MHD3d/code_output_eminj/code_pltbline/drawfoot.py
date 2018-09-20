import numpy as np
import matplotlib.pyplot as plt
import sys

#rawxp = np.genfromtxt("../../../../binfile/coord.xgc")
#rawyp = np.genfromtxt("../../../../binfile/coord.ygc")
rawxp = np.genfromtxt("../../../../3dbin/coord.xgc")
rawyp = np.genfromtxt("../../../../3dbin/coord.ygc")
XP, YP = np.meshgrid(rawxp,rawyp)
XPC, YPC = np.meshgrid(rawxp,rawyp)
print(rawxp.shape)
infos = np.genfromtxt("../../../../code_output/info",
                      dtype="string",
                      comments="!")
tsini=int(infos[0]); tsfin=int(infos[1])
nx = int(infos[3]) - int(infos[2]) + 1
ny = int(infos[5]) - int(infos[4]) + 1
nz = int(infos[7]) - int(infos[6]) + 1
zcut = 0
name = "Tw2dbin."
head = ("head","<i4"); tail = ("tail","<i4")
sizebin = nx*ny; chrsizebin = str(sizebin)
raw = ("val",chrsizebin+"f4"); dt = np.dtype([head,raw,tail])
array_foot = np.zeros([tsfin-tsini+1,2])

for i in range(tsini,tsfin+1):
    ts = i; cts = str(ts).zfill(3)
    data = np.genfromtxt("diff_L"+cts).reshape(nx,ny).T
#data = np.genfromtxt("diff_L018").reshape(ny+1,nx+1)
    print(data.shape)
    fd = open("../../"+name+cts,"r")
    chunk = np.fromfile(fd,dtype=dt,count=1)
    val = chunk[0]["val"].reshape((nx,ny),order="F")
    
    plt.clf()
#    plt.figure(figsize=(12,6))
    plt.pcolormesh(XP,YP,data,vmin=-1.0,vmax=1.0,cmap="seismic")
    plt.colorbar()
#    CR = plt.contour(XPC,YPC,val.T,[-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0],colors="k")
#    CR = plt.contour(XPC,YPC,val.T,[-1.0,-0.5,0.0],colors=["k","r","blue"],linewidth=0.1)
#    CR = plt.contour(XPC,YPC,val.T)
#    plt.clabel(CR)
#    CR.clabel(CR,inline='False')
    plt.xlim([-4.0,4.0])
    plt.ylim([-4.0,4.0])
    plt.axes().set_aspect("equal","datalim")
######################################################
    leng_crit = 0.1
    tmpix = 0
    tmpjy = 0
    for ix in range(0,nx):
        for jy in range(0,ny):
            if data[jy,ix] > leng_crit:
                if rawxp[ix] > rawxp[tmpix]:
                    if abs(rawxp[ix]) < 4.0:
                        if abs(rawyp[jy]) < 4.0:
                            tmpix = ix
                            tmpjy = jy
    array_foot[i-tsini,0] = tmpix
    array_foot[i-tsini,1] = tmpjy
######################################################
    plt.text(-4.0,5.0,"x:"+str(rawxp[tmpix]))
    plt.text(-4.0,4.6,"y:"+str(rawyp[tmpjy]))
    plt.text(-4.0,4.2,"L:"+str(data[tmpjy,tmpix]))
    plt.plot(rawxp[tmpix],rawyp[tmpjy],marker="+")
    plt.savefig("footdiff"+cts)
###
    np.savetxt("foot_edge.dat",array_foot)
