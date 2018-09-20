import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

rawxp = np.genfromtxt("../../binfile/coord.xgc")
rawyp = np.genfromtxt("../../binfile/coord.ygc")
rawzp = np.genfromtxt("../../binfile/coord.zgc")
rawtime1 = np.genfromtxt("../../../TIME_02")
time1 = rawtime1[:,2]
infos = np.loadtxt("../../code_output/info",comments="!")

tsini = int(infos[0])
tsfin = int(infos[1])
ts = 2
timeinit = 61.0
timefin = 100.0
time = np.append(timeinit,time1)
time = np.append(time,timefin)
cts = str(ts).zfill(3)
zcut = 0
zdrawcb = 190
ydrawcbmin = 256-200
ydrawcbmax = 256+200
xdrawcbmin = 256-200
xdrawcbmax = 256+200
#xcut = 35
xcut = 256
ycut = 256
nnx = 513
nny = 513
nnz = 257
vars = ["Bx","By","Bz","Vx","Vy","Vz","Phi"]
#vars = ["Vx","Vy","Vz","Bx","By","Bz","Ro"]
Bmax = 1.0
Bmin =-1.0
Vmax = 0.1
Vmin =-0.1

head = ("head","<i4")
tail = ("tail","<i4")
sizebin = nnx*nny*nnz
chrsizebin = str(sizebin)
raw = ("data",chrsizebin+"<f4")
dt = np.dtype([head,raw,tail])
n0xp = 0
n1xp = nnx+1
n0yp = 0
n1yp = nny+1
limiter_rho = np.ones([nnx,nny,nnz])*0.00000000000001
print(limiter_rho)

for var in vars:
    print("DRAW :",var)
    fi = open(var+"/check.txt","w")
    for ts in range(tsini,tsfin+1):
        print("---",ts)
        fi.write("ts:"+str(ts)+"\n")
        cts = str(ts).zfill(3)
        fd = open("../../binfile/"+var+"bin_3d."+cts,'r')
        chunk = np.fromfile(fd,dtype=dt,count=1)
        data = chunk[0]["data"].reshape((nnx,nny,nnz),order="F")
        fd.close()

#    print(data[:,:,0].shape)
#    print(rawxp.shape)
        XP,YP = np.meshgrid(rawxp,rawyp)
        YP_x,ZP_x = np.meshgrid(rawyp,rawzp) 
        XP_y, ZP_y = np.meshgrid(rawxp,rawzp)
    #    print(XP.shape)
        print(data[:,:,:].argmax())
        print(np.unravel_index(data.argmax(),data.shape))
    #   print(np.unravel_index(data.argmax(),data.shape)[1])
        tmp0 = np.unravel_index(data.argmax(),data.shape)[0]
        tmp1 = np.unravel_index(data.argmax(),data.shape)[1]
        tmp2 = np.unravel_index(data.argmax(),data.shape)[2]
        print("maxarr: ", tmp0, tmp1, tmp2)
        print("maxval: ", data[tmp0,tmp1,tmp2])
        fi.write("maxarr:"+str(tmp0)+","+str(tmp1)+","+str(tmp2)+"\n")
        fi.write("maxval:"+str(data[tmp0,tmp1,tmp2])+"\n")
        mintmp0 = np.unravel_index(data.argmin(),data.shape)[0]
        mintmp1 = np.unravel_index(data.argmin(),data.shape)[1]
        mintmp2 = np.unravel_index(data.argmin(),data.shape)[2]
        print("minarr: ", mintmp0, mintmp1, mintmp2)
        print("minval: ", data[mintmp0,mintmp1,mintmp2])
        fi.write("minarr:"+str(mintmp0)+","+str(mintmp1)+","+str(mintmp2)+"\n")
        fi.write("minval:"+str(data[mintmp0,mintmp1,mintmp2])+"\n")
        
#        plt.clf()
        plt.close("all")
        fig = plt.figure()
        ax = fig.add_subplot(111)
#        plt.pcolormesh(XP,YP,data[:,:,zcut].T)
#        plt.pcolormesh(YP_x,ZP_x,data[xcut,:,:].T)
        if var == "Ro":
#            plt.pcolormesh(YP_x,ZP_x,data[xcut,:,:].T)
            data[xcut,:,:] = np.maximum(data[xcut,:,:],limiter_rho[xcut,:,:])
            im = ax.pcolormesh(YP_x,ZP_x,data[xcut,:,:].T,norm=colors.LogNorm(vmin=np.abs(data[xcut,:,:].min()),vmax=data[xcut,:,:].max()))
        else:
#            plt.pcolormesh(YP_x,ZP_x,data[xcut,:,:].T)
            im = ax.pcolormesh(XP,YP,data[:,:,zcut].T,
                           vmin=np.min(data[xdrawcbmin:xdrawcbmax,
                                            ydrawcbmin:ydrawcbmax,
                                            zcut]),
                           vmax=np.max(data[xdrawcbmin:xdrawcbmax,
                                            ydrawcbmin:ydrawcbmax,
                                            zcut]),
                           cmap="bwr")
#        plt.pcolormesh(XP_y,ZP_y,data[:,ycut,:].T)
#            divider = make_axes_locatable(ax)
#            cax = divider.append_axes("right",size="5%",pad=0.05)
            plt.colorbar(im)
#
        ax.contour(XP,YP,data[:,:,zcut].T,
                   levels=[0.0],
                   linestyles="dashed"
                   )
        ax.set_title("color map "+var+" on xy = "+str(zcut))
        ax.set_xlabel("x")
        ax.set_ylabel("y")
#
#######     drawing bottom region
        ax.text(-4.0,-4.8,"t_A="+str(time[ts-tsini]))
        ax.set_xlim([-4.0,4.0])
        ax.set_ylim([-4.0,4.0])
        ax.set_aspect("equal")
        plt.savefig(var+"/bottom_"+var+str(zcut)+"yz"+cts+".png",dpi=300)
#######     drawing all region
        # ax.text(-75.0,-15.0,"t_A="+str(time[ts-tsini]))
#        ax.set_aspect("equal")
        # plt.savefig(var+"/all_"+var+str(xcut)+"yz"+cts+".png",dpi=300)
#######
#plt.show()

    fi.close()
