import numpy as np
import matplotlib.pyplot as plt
import sys
import glob as gb

rawxp = np.genfromtxt("../3dbin/coord.xgc")
rawyp = np.genfromtxt("../3dbin/coord.ygc")
rawzp = np.genfromtxt("../3dbin/coord.zgc")
strtime = gb.glob("../../TIME*")
rawtime = np.genfromtxt(strtime[0])
time = rawtime[:,2]
infos = np.loadtxt("../code_output/info",
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

timefin = get_info_restartfile("../../")
timeinit = get_info_restartfile(infos[8])
time = np.append(timeinit,time)
time = np.append(time,timefin)
tsini = int(infos[0]); tsfin = int(infos[1])
nnx = int(infos[3]) - int(infos[2]) + 1
nny = int(infos[5]) - int(infos[4]) + 1
nnz = int(infos[7]) - int(infos[6]) + 1
cts = str(tsini).zfill(3)
xcut = nnx/2; ycut = nny/2; zcut = 0
var = "By"

head = ("head","<i4"); tail = ("tail","<i4")
sizebin = nnx*nny*nnz; chrsizebin = str(sizebin)
raw = ("data",chrsizebin+"<f4"); dt = np.dtype([head,raw,tail])
sizebin2d = nnx*nny; chrsizebin2d = str(sizebin2d)
raw2d = ("data2d",chrsizebin2d+"<f4"); dt2d = np.dtype([head,raw2d,tail])

fd = open("../binfile/"+var+"bin_3d."+cts,'r')
chunk = np.fromfile(fd,dtype=dt,count=1)
data = chunk[0]["data"].reshape((nnx,nny,nnz),order="F")

#XP,YP = np.meshgrid(rawxp,rawyp)
XP, ZP = np.meshgrid(rawxp,rawzp)

dx = np.zeros(rawxp.shape[0])
dy = np.zeros(rawyp.shape[0])
dz = np.zeros(rawzp.shape[0])
#
for i in range(0,nnx):
    if i == 0:
        dx[i] = rawxp[i+1] - rawxp[i]
    elif i == nnx-1:
        dx[i] = rawxp[i] - rawxp[i-1]
    else:
        dx[i] = (rawxp[i+1] - rawxp[i-1])*0.5
#
for j in range(0,nny):
    if j == 0:
        dy[j] = rawyp[j+1] - rawyp[j]
    elif j == nny-1:
        dy[j] = rawyp[j] - rawyp[j-1]
    else:
        dy[j] = (rawyp[j+1] - rawyp[j-1])*0.5
#
for k in range(0,nnz):
    if k == 0:
        dz[k] = rawzp[k+1] - rawzp[k]
    elif k == nnz-1:
        dz[k] = rawzp[k] - rawzp[k-1]
    else:
        dz[k] = (rawzp[k+1] - rawzp[k-1])*0.5
#
segxz = np.zeros((rawxp.shape[0],rawzp.shape[0]))
segxy = np.zeros((rawxp.shape[0],rawyp.shape[0]))
print(dx.shape)
print(segxz.shape)

#zadj = 0.039035
zadj = 0.01
tmpfoot = np.genfromtxt("../r3dbline/file_output/picture/foot/seedall.txt")
print("check",tmpfoot.shape[0])
print(tmpfoot[0,:])
print(tmpfoot[1,:])
idx = np.abs(rawxp - tmpfoot[0,0]).argmin()
print(idx)
idx = np.abs(rawxp - tmpfoot[1,0]).argmin()
print(idx)
tmpfoot[:,3] = tmpfoot[:,3] + tsini

filenum = tsfin - tsini + 1
kappa = np.zeros(filenum); kappa_above = np.zeros(filenum)
Phi_tot = np.zeros(filenum)
tmp0 = np.zeros(filenum); tmp1 = np.zeros(filenum)
tst = 0

#### calculate Phi_total 
for ts_pt in range(tsini,tsfin+1):
    cts = str(ts_pt).zfill(3)
    fd = open("../3dbin/Bybin_3d_R."+cts,'r')
    chunk = np.fromfile(fd,dtype=dt,count=1)
    databy = chunk[0]["data"].reshape((nnx,nny,nnz),order="F")
    for i in range(0,nnx):
        for k in range(0,nnz):
            #
            segxz[i,k] = dx[i]*dz[k]
            Phi_tot[ts_pt-tsini] = Phi_tot[ts_pt-tsini] + segxz[i,k] * databy[i,ycut,k]

#### make array segxy
for i in range(0,nnx):
    for j in range(0,nny):
        segxy[i,j] = dx[i]*dy[j]

# ######### debug cal by Phi_tot_on bz
Phi_tot_test = np.zeros(filenum)
zcut = 0
for ts_pt in range(tsini,tsfin+1):
    cts = str(ts_pt).zfill(3)
    fd = open("../3dbin/Bzbin_3d_R."+cts,'r')
    chunk = np.fromfile(fd,dtype=dt,count=1)
    databz = chunk[0]["data"].reshape((nnx,nny,nnz),order="F")
    for i in range(0,nnx):
        for j in range(0,nny):
            #
            Phi_tot_test[ts_pt-tsini] = Phi_tot_test[ts_pt-tsini] + segxy[i,j] * np.abs(databz[i,j,zcut])*0.5
print(Phi_tot_test)
#sys.exit()

#### calculate \int dPhi_rec * Twist
for tt in range(0,tmpfoot.shape[0]):
    ix = np.abs(rawxp - tmpfoot[tt,0]).argmin()
    jy = np.abs(rawyp - tmpfoot[tt,1]).argmin()
    kz = np.abs(rawzp - tmpfoot[tt,2]).argmin()
    if(int(tmpfoot[tt,3])!=tst):
        tst = int(tmpfoot[tt,3])
        cts = str(int(tst)).zfill(3)
        #
        fdbz = open("../3dbin/Bzbin_3d_R."+cts,'r')
        chunkbz = np.fromfile(fdbz,dtype=dt,count=1)
        databz = chunkbz[0]["data"].reshape((nnx,nny,nnz),order="F")
        #
        fdtw = open("../r3dbline/file_output/Tw2dbin."+cts,'r')
        chunktw = np.fromfile(fdtw,dtype=dt2d,count=1)
        datatw = chunktw[0]["data2d"].reshape((nnx,nny),order="F")
        
#    print(tst,ix,jy,kz)
    kappa_above[tst-tsini] = kappa_above[tst-tsini] + datatw[ix,jy]*segxy[ix,jy]*abs(databz[ix,jy,0])*0.5
#    kappa_above[tst] = kappa_above[tst] + datatw[ix,jy]*segxy[ix,jy]*abs(databz[ix,jy,0])*0.5

kappa = kappa_above/Phi_tot
kappa = np.abs(kappa)
#print(kappa_above)
#print(Phi_tot)
#print(kappa)
kappa[np.isnan(kappa)] = 0
#print(kappa)

print(time)
print(time.shape)
print(kappa)
print(kappa.shape)

### make figure
### plot for kappa
plt.clf()
fig = plt.figure()
fig.subplots_adjust(left=0.18,right=0.96)
ax = fig.add_subplot(111)
#plt.plot(time[tsini:tsfin],kappa[:-1])
#plt.plot([time[tsini-tsini],time[tsfin-2]],[1.0/(4.0*np.pi),1.0/(4.0*np.pi)],linestyle="--",dashes=(12,6))
#plt.plot([time[tsini],time[tsfin-2]],[1.0/8.0,1.0/8.0],linestyle="--",dashes=(12,6,4,6))
#plt.plot([time[tsini],time[tsfin-2]],[7.0/40.0,7.0/40.0],linestyle="--",dashes=(12,6,4,6,4,6))
plt.plot(time[:],kappa[:])
plt.plot([time[tsini-tsini],time[tsfin-tsini]],[1.0/(4.0*np.pi),1.0/(4.0*np.pi)],linestyle="--",dashes=(12,6))
plt.plot([time[tsini-tsini],time[tsfin-tsini]],[1.0/8.0,1.0/8.0],linestyle="--",dashes=(12,6,4,6))
plt.plot([time[tsini-tsini],time[tsfin-tsini]],[7.0/40.0,7.0/40.0],linestyle="--",dashes=(12,6,4,6,4,6))
#plt.plot(abs(Phi_tot))
#plt.plot(abs(kappa_above))
plt.xlabel("time")
plt.ylabel("$\kappa$",fontsize=20)
#plt.xlim([0,30])
#plt.show()
plt.savefig("kappa.png")

### plot for Phi_total
plt.clf()
fig = plt.figure()
fig.subplots_adjust(left=0.18,right=0.96)
ax = fig.add_subplot(111)
#plt.plot(time[tsini:tsfin],Phi_tot[:-1])
plt.plot(time[:],Phi_tot[:])
plt.xlabel("time")
plt.ylabel("$\Phi _{tot}$",fontsize=20)
#plt.show()
plt.savefig("Phi_tot.png")

### plot for \int T dphi
plt.clf()
fig = plt.figure()
fig.subplots_adjust(left=0.18,right=0.96)
ax = fig.add_subplot(111)
#plt.plot(time[tsini:tsfin],kappa_above[:-1])
plt.plot(time[:],kappa_above[:])
plt.xlabel("time")
#plt.ylabel("$\kappa _{numerator}$",fontsize=20)
plt.ylabel("$\int_{rec} T d\phi$",fontsize=20)
#plt.show()
plt.savefig("kappa_above.png")

np.savetxt("kappa_t"+str(timeinit)+".dat",
           [time[:],kappa[:]])
