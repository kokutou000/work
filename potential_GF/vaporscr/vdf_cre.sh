#!/bin/sh
FILENAME="val.vdf"
COORDDIR="../"
DATADIR="./"
XPARA=514
YPARA=514
ZPARA=257
vdfcreate -dimension ${XPARA}x${YPARA}x${ZPARA} -xcoords ${COORDDIR}coord.xgc -ycoords ${COORDDIR}coord.ygc -zcoords ${COORDDIR}coord.zgc -numts 1 -vars3d bx:by:bz ${FILENAME} -gridtype stretched

raw2vdf -ts 0 -varname bx ${FILENAME} ${DATADIR}Bx4bin3d
raw2vdf -ts 0 -varname by ${FILENAME} ${DATADIR}By4bin3d
raw2vdf -ts 0 -varname bz ${FILENAME} ${DATADIR}Bz4bin3d
