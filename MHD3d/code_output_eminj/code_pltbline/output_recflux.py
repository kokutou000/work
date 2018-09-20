import numpy as np
import matplotlib.pyplot as plt
import glob as gb
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__),
                             "../../../../code_output"))
infos = np.loadtxt("../../../../code_output/info",
                   dtype="string",comments="!")
strtime = gb.glob("../../../../../TIME*")
rawtime = np.genfromtxt(strtime[0])
import get_info_restart as gir
timeinit=gir.get_info_restartfile("../../../"+infos[8],"LITTLE_ENDIAN")
timefin=gir.get_info_restartfile("../../../../../","BIG_ENDIAN")
time = rawtime[:,2]
time = np.append(timeinit,time)
time = np.append(time,timefin)
infos = np.loadtxt("../../../../code_output/info",
                   dtype="string",comments="!")
endt_inj = float(infos[10])

rawrecfluxdata = np.genfromtxt("recflux.txt")
print(rawrecfluxdata.shape)
print(time.shape)

plt.plot(time,rawrecfluxdata[:,1],marker="+")
plt.axvline(x=endt_inj,linestyle="dashed",color="black")
plt.xlabel("time")
plt.ylabel("reconnected flux")
#plt.show()
plt.savefig("fig_recflux.png")
