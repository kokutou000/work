import numpy as np
import matplotlib.pyplot as plt
import sys

nnx = 514
nny = 514
tmpdbzdx = np.genfromtxt("dbzdx.dat")
tmpdbzdy = np.genfromtxt("dbzdy.dat")
tmpbz0 = np.genfromtxt("../Bz_init.dat")
tmpdbzdytanh = np.genfromtxt("dbzdytanhbz.dat")
tmpdbzdycal = np.genfromtxt("dbzdycal.dat")
tmpvxtmp = np.genfromtxt("vxtmp.dat")
tmpvytmp = np.genfromtxt("vytmp.dat")

print(tmpdbzdx.shape)
print(nnx*nny)
dbzdx = tmpdbzdx[:,2].reshape(nnx,nny)
dbzdy = tmpdbzdy[:,2].reshape(nnx,nny)
bz0 = tmpbz0[:,2].reshape(nnx,nny)
dbzdytanh = tmpdbzdytanh[:,2].reshape(nnx,nny)
dbzdycal = tmpdbzdycal[:,1].reshape(nnx,nny)
vxtmp = tmpvxtmp.reshape(nnx,nny)
vytmp = tmpvytmp.reshape(nnx,nny)
print(dbzdx.shape)

xc = np.genfromtxt("../coord.xgc")
yc = np.genfromtxt("../coord.ygc")
XC, YC = np.meshgrid(xc,yc)

vvnorm = np.max(np.sqrt(vxtmp**2+vytmp**2))
print(vvnorm)

plt.close("all")
plt.pcolormesh(XC,YC,(vxtmp/vvnorm).T,cmap="jet")
#plt.pcolormesh(XC,YC,vxtmp.T,cmap="jet")
plt.colorbar()
plt.xlim([-5.0,5.0])
plt.ylim([-5.0,5.0])
plt.title("vx z=0 normalized by max|v|")
plt.savefig("vxmap_modified.png")

plt.close("all")
plt.pcolormesh(XC,YC,(vytmp/vvnorm).T,cmap="jet")
#plt.pcolormesh(XC,YC,vytmp.T,cmap="jet")
plt.colorbar()
plt.xlim([-5.0,5.0])
plt.ylim([-5.0,5.0])
plt.title("vy z=0 normalized by max|v|")
plt.savefig("vymap_modified.png")

sys.exit()

plt.close("all")
plt.plot(yc,bz0[257,:],marker="+")
plt.plot(yc,dbzdy[257,:]/(np.max(dbzdy[257,:])),marker="+")
plt.text(0.0,1.1,"max(dbzdy):"+str(np.max(dbzdy[257,:])))
plt.xlim([0.0,5.0])
plt.ylim([-0.5,1.0])
plt.savefig("bz_dbzdy__vs__y.png")

plt.close("all")
plt.plot(yc,bz0[257,:],marker="+")
plt.plot(yc,dbzdy[257,:]/(np.max(dbzdy[257,:])),marker="+")
plt.text(0.0,1.1,"max(dbzdy):"+str(np.max(dbzdy[257,:])))
plt.plot(yc,dbzdytanh[257,:]/(np.max(dbzdytanh[257,:])),marker="^")
plt.text(0.0,1.05,"max(dbzdytanh):"+str(np.max(dbzdytanh[257,:])))
plt.xlim([0.0,5.0])
plt.ylim([-0.5,1.0])
plt.savefig("dbzdytanh__vs__y.png")


plt.close("all")
plt.plot(yc,bz0[257,:],marker="+")
plt.plot(yc,dbzdy[257,:]/(np.max(dbzdy[257,:])),marker="+")
plt.text(0.0,1.1,"max(dbzdy):"+str(np.max(dbzdy[257,:])))
plt.plot(yc,dbzdytanh[257,:]/(np.max(dbzdytanh[257,:])),marker="^")
plt.text(0.0,1.05,"max(dbzdytanh):"+str(np.max(dbzdytanh[257,:])))
plt.plot(yc,vxtmp[257,:]/vvnorm,marker="+")
plt.xlim([0.0,5.0])
plt.ylim([-0.5,1.0])
plt.savefig("vxs__vs__y.png")


### make plot of normalizing function ###
tmpnewfunc = np.genfromtxt("newfunction.dat")
plt.close("all")
plt.plot(tmpdbzdycal[:,0],tmpdbzdycal[:,1],
         label="dbzdy_vs_bz",linestyle="None",
         marker=".")
plt.plot(tmpnewfunc[:,0],tmpnewfunc[:,1],
         label="new function",linestyle="None",
         marker=".")
plt.xlim([0.0,1.0])
plt.xlabel("Bz")
plt.ylabel("f(Bz) or dBzdy(Bz)")
plt.legend(loc="upper right")
plt.savefig("plot_functions.png")
