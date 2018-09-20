#!/bin/sh
FILENAME="val.vdf"
COORDDIR="../3dbin/"
DATADIR="../3dbin/"
DATA2DIR="../r3dbline/file_output/"
IFS=$'\n'
file=(`cat ../code_output/info`)
XPARA=$((file[3]-file[2]+1))
YPARA=$((file[5]-file[4]+1))
ZPARA=$((file[7]-file[6]+1))
numb_ts=1
num_fints=${file[1]}
echo ${num_fints}
###
vdfcreate -dimension ${XPARA}x${YPARA}x${ZPARA} -xcoords ${COORDDIR}coord.xgc -ycoords ${COORDDIR}coord.ygc -zcoords ${COORDDIR}coord.zgc -numts ${numb_ts} -vars3d bx:by:bz:vx:vy:vz:jx:jy:jz:jj ${FILENAME} -vars2dxy bz2d:tw2d:al2d -gridtype stretched
raw2vdf -ts 0 -varname bx ${FILENAME} ${DATADIR}Bxbin_3d_R.$(printf "%03d" ${num_fints})
raw2vdf -ts 0 -varname by ${FILENAME} ${DATADIR}Bybin_3d_R.$(printf "%03d" ${num_fints})
raw2vdf -ts 0 -varname bz ${FILENAME} ${DATADIR}Bzbin_3d_R.$(printf "%03d" ${num_fints})
raw2vdf -ts 0 -varname vx ${FILENAME} ${DATADIR}Vxbin_3d_R.$(printf "%03d" ${num_fints})
raw2vdf -ts 0 -varname vy ${FILENAME} ${DATADIR}Vybin_3d_R.$(printf "%03d" ${num_fints})
raw2vdf -ts 0 -varname vz ${FILENAME} ${DATADIR}Vzbin_3d_R.$(printf "%03d" ${num_fints})
raw2vdf -ts 0 -varname jx ${FILENAME} ${DATA2DIR}Jx3dbin.$(printf "%03d" ${num_fints})
raw2vdf -ts 0 -varname jy ${FILENAME} ${DATA2DIR}Jy3dbin.$(printf "%03d" ${num_fints})
raw2vdf -ts 0 -varname jz ${FILENAME} ${DATA2DIR}Jz3dbin.$(printf "%03d" ${num_fints})
raw2vdf -ts 0 -varname jj ${FILENAME} ${DATA2DIR}JJ3dbin.$(printf "%03d" ${num_fints})
raw2vdf -ts 0 -varname bz2d ${FILENAME} ${DATADIR}Bz2d_R.$(printf "%03d" ${num_fints})
raw2vdf -ts 0 -varname tw2d ${FILENAME} ${DATA2DIR}Tw2dbin.$(printf "%03d" ${num_fints})
raw2vdf -ts 0 -varname al2d ${FILENAME} ${DATA2DIR}alp2dbin.$(printf "%03d" ${num_fints})
exit
