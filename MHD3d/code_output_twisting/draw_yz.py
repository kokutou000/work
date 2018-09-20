import numpy as np
import matplotlib.pyplot as plt
import sys

rawxp = np.genfromtxt("../../binfile/coord.xgc")
rawyp = np.genfromtxt("../../binfile/coord.ygc")
rawzp = np.genfromtxt("../../binfile/coord.zgc")
rawtime = np.genfromtxt("../../../TIME_02")
infos = np.loadtxt("../../code_output/info",comments="!")

tsini = int(infos[0])
tsfin = int(infos[1])
ts = 2
cts = str(ts).zfill(3)
zcut = 0
#xcut = 35
xcut = 256
ycut = 256
timeinit = 200.0
nnx = 513
nny = 513
nnz = 257
nnxlow = 88
nnxupp = 424
nnylow = 88
nnyupp = 424
nnzlow = 0
nnzupp = 163
vars = ["Bx","By","Bz","Vx","Vy","Vz","Phi"]
#vars = ["Vx","Vy","Vz"]
Bmax = 1.0
Bmin =-1.0
Vmax = 0.1
Vmin =-0.1
timearr = np.array([timeinit])
timearr = np.hstack((timearr,rawtime[:,2]))
#timearr = np.append(rawtime[:,2],rawtime2[:,2])
timearr = np.hstack((timearr,np.array([999.9])))
timearr = np.hstack((np.array(timeinit),timearr))
n_yz = np.genfromtxt("../../../../decay_index/n_yz.dat")
n_yz = n_yz.reshape(nny+1,nnz)
print(n_yz.shape)
n_yz = n_yz[1:,:]
print(n_yz.shape)

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

for var in vars:
    print("DRAW :",var)
    fi = open(var+"/check.txt","w")
    for ts in range(tsini,tsfin+1):
        print("---",ts)
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
        print(data[tmp0,tmp1,tmp2])
        fi.write("maxarr:"+str(tmp0)+","+str(tmp1)+","+str(tmp2)+"\n")
        fi.write("maxval:"+str(data[tmp0,tmp1,tmp2])+"\n")
        mintmp0 = np.unravel_index(data.argmin(),data.shape)[0]
        mintmp1 = np.unravel_index(data.argmin(),data.shape)[1]
        mintmp2 = np.unravel_index(data.argmin(),data.shape)[2]
        fi.write("minarr:"+str(mintmp0)+","+str(mintmp1)+","+str(mintmp2)+"\n")
        fi.write("minval:"+str(data[mintmp0,mintmp1,mintmp2])+"\n")
        
        plt.clf()
        fig = plt.figure()
        ax = fig.add_subplot(111)
#        plt.pcolormesh(XP,YP,data[:,:,zcut].T)
#        plt.pcolormesh(YP_x,ZP_x,data[xcut,:,:].T)
### plot bottom region
        im = ax.pcolormesh(YP_x,ZP_x,data[xcut,:,:].T,
                           vmax=np.max(data[xcut,
                                            nnylow:nnyupp,
                                            nnzlow:nnzupp]),
                           vmin=np.min(data[xcut,
                                            nnylow:nnyupp,
                                            nnzlow:nnzupp]),
                           cmap="bwr"
                           )
        #    plt.pcolormesh(XP_y,ZP_y,data[:,ycut,:].T)
        plt.colorbar(im,orientation="horizontal")
        ax.contour(YP_x,ZP_x,data[xcut,:,:].T,
                   levels=[0.0],
                   linestyles="dashed"
                   )
#### plot contour of decay index
        ax.contour(YP_x,ZP_x,n_yz.T,
                  levels=[0.5,1.0,1.5],
                  linestyles="dashed"
                  )
        ax.set_title(var+" on yz plame x="+str(xcut))
        ax.text(-4.0,-0.8,"t = "+str(timearr[ts-tsini])+" t_A")
        ax.set_xlim([-2.0,2.0])
        ax.set_ylim([0.0,2.0])
        ax.set_aspect("equal")
        plt.savefig(var+"/bottom"+var+"yz"+cts+".png",dpi=300)
### plot all region
#         im = ax.pcolormesh(YP_x,ZP_x,data[xcut,:,:].T,cmap="bwr")
#         plt.colorbar(im)
#         plt.contour(YP_x,ZP_x,data[xcut,:,:].T,
#                     levels=[0.0],
#                     linestyles="dashed"
#                     )
#         ax.set_title(var+" on yz plame x="+str(xcut))
#         ax.text(-75.0,-15.0,"t = "+str(timearr[ts-tsini])+" t_A")
# #        ax.set_xlim([-3.0,3.0])
# #        ax.set_ylim([0.0,6.0])
#         plt.savefig(var+"/all"+var+"yz"+cts+".png",dpi=300)
#####################
        plt.close()
