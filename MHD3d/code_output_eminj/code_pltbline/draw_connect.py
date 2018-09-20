import numpy as np
import matplotlib.pyplot as plt
import sys

infos = np.genfromtxt("../../../../code_output/info",
                      dtype="string",
                      comments="!")
tsini = int(infos[0]); tsfin = int(infos[1])
nx = int(infos[3]) - int(infos[2]) + 1
ny = int(infos[5]) - int(infos[4]) + 1

for ts in range(tsini,tsfin+1):
    cts = str(ts).zfill(3)
    tmp = np.genfromtxt("../../Byzp"+cts)
    
    originxy = tmp[0,:]
    connectyz = tmp[1,:]
    for i in range(1,tmp.shape[0]/2):
        originxy = np.vstack((originxy,tmp[i*2,:]))
        connectyz = np.vstack((connectyz,tmp[i*2+1,:]))
        
    print("orig.shape",originxy.shape)
    print("conn.shape",connectyz.shape)
    
    plt.scatter(connectyz[:,1],connectyz[:,2],
                marker=".",
                s=5)
    plt.savefig("yz_posi"+cts+"all.png")
    plt.close("all")
#    plt.show()
    
    diff = np.genfromtxt("../foot/diff_L"+cts)
    print("diff.shape",diff.shape)
    rawxp = np.genfromtxt("../../../../3dbin/coord.xgc")
    rawyp = np.genfromtxt("../../../../3dbin/coord.ygc")
    diff = diff.reshape(nx,ny)
    print("diff.shape",diff.shape)
    # plt.pcolormesh(diff.T)
    # plt.colorbar()
    # plt.show()
#
    for ix in range(0,rawxp.shape[0]):
        if round(originxy[0,0],5) == round(rawxp[ix],5):
            tmpix = np.array(ix)
#
    for iori in range(1,originxy.shape[0]):
        for ix in range(0,rawxp.shape[0]):
            if round(originxy[iori,0],5) == round(rawxp[ix],5):
                tmpix = np.append(tmpix,ix)
    print("tmpix",tmpix)
    print("tmpix.shape",tmpix.shape)
#
    for jy in range(0,rawyp.shape[0]):
        if round(originxy[0,1],5) == round(rawyp[jy],5):
            tmpjy = np.array(jy)
#
    for jori in range(1,originxy.shape[0]):
        for jy in range(0,rawyp.shape[0]):
            if round(originxy[jori,1],5) == round(rawyp[jy],5):
                tmpjy = np.append(tmpjy,jy)
    print("tmpjy",tmpjy)
    print("tmpjy.shape",tmpjy.shape)

#    plt.scatter(rawxp[tmpix],rawyp[tmpjy])
#    plt.scatter(tmpix,tmpjy)
#    plt.show()

    if diff[tmpix[0],tmpjy[0]] < 0.4:
        arr_yzrec = 0.0
    else:
        arr_yzrec = 1.0

    for na in range(1,tmpix.shape[0]):
        if diff[tmpix[na],tmpjy[na]] < 0.4:
            arr_yzrec = np.append(arr_yzrec,0)
        else:
            arr_yzrec = np.append(arr_yzrec,1)
    print("arr_yzrec",arr_yzrec)
    print("arr_yzrec.shape",arr_yzrec.shape)
    
    arr_draw = np.vstack((connectyz[:,1],connectyz[:,2]))
    arr_draw = np.vstack((arr_draw,arr_yzrec))
    print("arr_draw.shape",arr_draw.shape)
    print(arr_draw)
#
    sumx = 0.0
    sumy = 0.0
    count = 0
    for na in range(0,tmpix.shape[0]):
        if arr_draw[2,na] == 1.0:
            sumx = sumx + arr_draw[0,na]
            sumy = sumy + arr_draw[1,na]
            count = count + 1
    if count != 0:
        avex = sumx/count
        avey = sumy/count
    else:
        avex = 0.0
        avey = 0.0
    print(avex,avey,sumx,sumy,count)
#
    if ts == tsini:
        zevol = 0.0
    else:
        zevol = np.append(zevol,avey)
    
#    plt.scatter(connectyz[:,1],connectyz[:,2],
    plt.scatter(arr_draw[0,:],arr_draw[1,:],
                c=arr_draw[2,:],
                marker=".",
                s=5)
    plt.plot(avex,avey,
             marker="+",
             color="r")
    plt.ylim([0.0,15.0])
    plt.text(-10.0,-3.5,str(avex)+","+str(avey))
#    plt.colorbar()
    plt.savefig("FR_posi_yz"+cts+".png")
    plt.close("all")
#    plt.show()
    np.savetxt("outFRposi"+cts,arr_draw)
### make figure of histgram of FR point for z-height
    arr_histz = 0.0
    for na in range(0,tmpix.shape[0]):
        if arr_draw[2,na] == 1.0:
            arr_histz = np.append(arr_histz,arr_draw[1,na])

    plt.hist(arr_histz,bins=100,range=(0.0,10.0))
    plt.axvline(x=avey,linestyle="dashed",color="red")
    plt.xlim([0.0,15.0])
    plt.savefig("hist_FRz"+cts+".png")
    plt.close("all")

np.savetxt("evol_FRz.dat",zevol)
