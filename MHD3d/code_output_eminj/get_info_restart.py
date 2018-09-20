import numpy as np
import sys

def get_info_restartfile(path_refile,cendian):
#-----    endian_ck = 0 --> big_endian
#-----    endian_ck = 1 --> little_endian
    if cendian == "LITTLE_ENDIAN":
        sendi="<"
    elif cendian == "BIG_ENDIAN":
        sendi=">"
    else:
        "chr_endian is bad value!!!"
        sys.exit()
    headRe = ("head",sendi+"i4"); iwriteRe = ("iwrite",sendi+"i4")
    nloopRe = ("nloop",sendi+"i4"); atimeRe = ("atime",sendi+"f8")
    dtstepRe = ("dtstep",sendi+"f8"); tailRe = ("tail",sendi+"i4")
    nxre=66; nyre=66; nzre=257
    sizebinRe = nxre*nyre*nzre; chrsizebinRe = str(sizebinRe)
    rawRe = ("datare",chrsizebinRe+sendi+"f4")
    dtRe = np.dtype([headRe,iwriteRe,nloopRe,
                     atimeRe,dtstepRe,rawRe,tailRe])
    #
    fdre = open(path_refile+"RESTART.0000",'r')
    chunkre = np.fromfile(fdre,dtype=dtRe,count=1)
    datatime = chunkre[0]["atime"]
    fdre.close()
    return datatime
