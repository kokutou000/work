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
num_inits=${file[0]}
num_fints=${file[1]}
numb_ts=$((file[1]-file[0]+1))
### create vdf dummy file ###
vdfcreate -dimension ${XPARA}x${YPARA}x${ZPARA} -xcoords ${COORDDIR}coord.xgc -ycoords ${COORDDIR}coord.ygc -zcoords ${COORDDIR}coord.zgc -numts ${numb_ts} -vars3d bx:by:bz:vx:vy:vz:jx:jy:jz:jj ${FILENAME} -vars2dxy bz2d:tw2d:al2d -gridtype stretched
### reading binary data in do-loop ###
for j in `seq $num_inits $num_fints`
do
i=$((j - num_inits))
printf "ts vapor:i = $i \n"
printf "ts file:j = $j \n"
printf "++++++++++\n"
printf "%03d\n" $j
printf ${DATADIR}"Bxbin_3d."$(printf "%03d" $j)"\n"
printf "++++++++++\n"
raw2vdf -ts $i -varname bx ${FILENAME} ${DATADIR}Bxbin_3d_R.$(printf "%03d" $j)
raw2vdf -ts $i -varname by ${FILENAME} ${DATADIR}Bybin_3d_R.$(printf "%03d" $j)
raw2vdf -ts $i -varname bz ${FILENAME} ${DATADIR}Bzbin_3d_R.$(printf "%03d" $j)
raw2vdf -ts $i -varname vx ${FILENAME} ${DATADIR}Vxbin_3d_R.$(printf "%03d" $j)
raw2vdf -ts $i -varname vy ${FILENAME} ${DATADIR}Vybin_3d_R.$(printf "%03d" $j)
raw2vdf -ts $i -varname vz ${FILENAME} ${DATADIR}Vzbin_3d_R.$(printf "%03d" $j)
raw2vdf -ts $i -varname jx ${FILENAME} ${DATA2DIR}Jx3dbin.$(printf "%03d" $j)
raw2vdf -ts $i -varname jy ${FILENAME} ${DATA2DIR}Jy3dbin.$(printf "%03d" $j)
raw2vdf -ts $i -varname jz ${FILENAME} ${DATA2DIR}Jz3dbin.$(printf "%03d" $j)
raw2vdf -ts $i -varname jj ${FILENAME} ${DATA2DIR}JJ3dbin.$(printf "%03d" $j)
raw2vdf -ts $i -varname bz2d ${FILENAME} ${DATADIR}Bz2d_R.$(printf "%03d" $j)
raw2vdf -ts $i -varname tw2d ${FILENAME} ${DATA2DIR}Tw2dbin.$(printf "%03d" $j)
done
printf "\n+++++++++++\n"
printf "\nfinish create vdf file!\n"