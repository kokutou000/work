import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob as gb

rawxp = np.genfromtxt("../../binfile/coord.xgc")
rawyp = np.genfromtxt("../../binfile/coord.ygc")
rawzp = np.genfromtxt("../../binfile/coord.zgc")
strtime = gb.glob("../../../TIME*")
rawtime = np.genfromtxt(strtime[0])
time = rawtime[:,2]
infos = np.loadtxt("../../code_output/info",
                   dtype="string",
                   comments="!")
def get_info_restartfile(path_refile):
    headRe = ("head","<i4"); iwriteRe = ("iwrite","<i4")
    nloopRe = ("nloop","<i4"); atimeRe = ("atime","<f8")
    dtstepRe = ("dtstep","<f8"); tailRe = ("tail","<i4")
    nxre=66; nyre=66; nzre=257
    sizebinRe = nxre*nyre*nzre; chrsizebinRe = str(sizebinRe)
    rawRe = ("datare",chrsizebinRe+"<f4")
#    dtRe = np.dtype([headRe,rawRe,tailRe])
    dtRe = np.dtype([headRe,iwriteRe,nloopRe,
                     atimeRe,dtstepRe,rawRe,tailRe])
    #
    fdre = open(path_refile+"RESTART.0001",'r')
    chunkre = np.fromfile(fdre,dtype=dtRe,count=1)
    datatime = chunkre[0]["atime"]
    fdre.close()
    return datatime

timefin = get_info_restartfile("../../../")
timeinit = get_info_restartfile("../"+infos[8])
tsini = int(infos[0])
tsfin = int(infos[1])
time = np.append(timeinit,time)
time = np.append(time,timefin)
xdrawl = int(infos[2]); xdrawu = int(infos[3])
ydrawl = int(infos[4]); ydrawu = int(infos[5])
zdrawl = int(infos[6]); zdrawu = int(infos[7])
xcut = 256; ycut = 256; zcut = 0
nnx = 513; nny = 513; nnz = 257
vars = ["Bx", "Vz", "Phi"]
#vars = ["Bx","By","Bz","Vx","Vy","Vz","Phi"]
#vars = ["Vx","Vy","Vz","Bx","By","Bz","Ro"]

### obtain and reform the data of decay index
n_yz = np.genfromtxt("../../../../decay_index/n_yz.dat")
n_yz = n_yz.reshape(nny+1,nnz)
n_yz = n_yz[1:,:]
###

head = ("head","<i4"); tail = ("tail","<i4")
sizebin = nnx*nny*nnz; chrsizebin = str(sizebin)
raw = ("data",chrsizebin+"<f4"); dt = np.dtype([head,raw,tail])
limiter_rho = np.ones([nnx,nny,nnz])*0.00000000000001

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

        XP,YP = np.meshgrid(rawxp,rawyp)
        YP_x,ZP_x = np.meshgrid(rawyp,rawzp) 
        XP_y, ZP_y = np.meshgrid(rawxp,rawzp)
        print(data[:,:,:].argmax())
        print(np.unravel_index(data.argmax(),data.shape))
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
        
        plt.close("all")
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if var == "Ro":
            data[xcut,:,:] = np.maximum(data[xcut,:,:],limiter_rho[xcut,:,:])
            im = ax.pcolormesh(YP_x,ZP_x,data[xcut,:,:].T,norm=colors.LogNorm(vmin=np.abs(data[xcut,:,:].min()),vmax=data[xcut,:,:].max()))
        else:
            im = ax.pcolormesh(YP_x,ZP_x,data[xcut,:,:].T,
                               vmin=np.min(data[xcut,
                                                ydrawl:ydrawu,
                                                zdrawl:zdrawu]),
                               vmax=np.max(data[xcut,
                                                ydrawl:ydrawu,
                                                zdrawl:zdrawu]),
                               cmap="bwr")
#            divider = make_axes_locatable(ax)
#            cax = divider.append_axes("right",size="5%",pad=0.05)
            plt.colorbar(im,orientation="horizontal")
#
        ax.contour(YP_x,ZP_x,data[xcut,:,:].T,
                   levels=[0.0],
                   linestyles="dashed"
                   )
### plot contour of decay index
        ax.contour(YP_x,ZP_x,n_yz.T,
                   levels=[0.5,1.0,1.5],
                   linestyles="dashed"
                   )
        ax.set_title("color map "+var+" on yz: x = "+str(xcut))
        ax.set_xlabel("y")
        ax.set_ylabel("z")

#######     drawing limited region
        ax.text(-2.0,-0.3,"t_A="+str(time[ts-tsini]))
        ax.set_xlim([-2.0,2.0])
        ax.set_ylim([0.0,2.0])
        ax.set_aspect("equal")
        plt.savefig(var+"/bottom_"+var+str(xcut)+"yz"+cts+".png",dpi=300)
#######     drawing all region
        # ax.text(-75.0,-15.0,"t_A="+str(time[ts-tsini]))
#        ax.set_aspect("equal")
        # plt.savefig(var+"/all_"+var+str(xcut)+"yz"+cts+".png",dpi=300)
#######
#plt.show()

    fi.close()
