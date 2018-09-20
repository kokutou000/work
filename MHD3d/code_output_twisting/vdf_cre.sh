#!/bin/sh
FILENAME="val.vdf"
COORDDIR="../3dbin/"
DATADIR="../3dbin/"
DATA2DIR="../r3dbline/file_output/"
XPARA=337
YPARA=337
ZPARA=164
numb_ts=41
###
vdfcreate -dimension ${XPARA}x${YPARA}x${ZPARA} -xcoords ${COORDDIR}coord.xgc -ycoords ${COORDDIR}coord.ygc -zcoords ${COORDDIR}coord.zgc -numts ${numb_ts} -vars3d bx:by:bz:vx:vy:vz:jx:jy:jz:jj ${FILENAME} -vars2dxy bz2d:tw2d:al2d -gridtype stretched
num=("001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" "015" "016" "017" "018" "019" "020" "021" "022" "023" "024" "025" "026" "027" "028" "029" "030" "031" "032" "033" "034" "035" "036" "037" "038" "039" "040" "041" "042" "043" "044" "045" "046" "047" "048" "049" "050" "051" "052" "053" "054" "055" "056" "057" "058" "059" "060" "061" "062" "063")
for i in {18..59}
do
printf "++++++++++\n"
printf ${num[$i]}"\n"
printf ${DATADIR}"Bxbin_3d."${num[$i]}"\n"
printf "++++++++++\n"
raw2vdf -ts $i -varname bx ${FILENAME} ${DATADIR}Bxbin_3d_R.${num[$i]}
raw2vdf -ts $i -varname by ${FILENAME} ${DATADIR}Bybin_3d_R.${num[$i]}
raw2vdf -ts $i -varname bz ${FILENAME} ${DATADIR}Bzbin_3d_R.${num[$i]}
raw2vdf -ts $i -varname vx ${FILENAME} ${DATADIR}Vxbin_3d_R.${num[$i]}
raw2vdf -ts $i -varname vy ${FILENAME} ${DATADIR}Vybin_3d_R.${num[$i]}
raw2vdf -ts $i -varname vz ${FILENAME} ${DATADIR}Vzbin_3d_R.${num[$i]}
raw2vdf -ts $i -varname jx ${FILENAME} ${DATA2DIR}Jx3dbin.${num[$i]}
raw2vdf -ts $i -varname jy ${FILENAME} ${DATA2DIR}Jy3dbin.${num[$i]}
raw2vdf -ts $i -varname jz ${FILENAME} ${DATA2DIR}Jz3dbin.${num[$i]}
raw2vdf -ts $i -varname jj ${FILENAME} ${DATA2DIR}JJ3dbin.${num[$i]}
raw2vdf -ts $i -varname bz2d ${FILENAME} ${DATADIR}Bz2d_R.${num[$i]}
raw2vdf -ts $i -varname tw2d ${FILENAME} ${DATA2DIR}Tw2dbin.${num[$i]}
raw2vdf -ts $i -varname al2d ${FILENAME} ${DATA2DIR}alp2dbin.${num[$i]}
done
printf "\n+++++++++++\n"
printf "\nfinish create vdf file!\n"
#raw2vdf -ts 0 -varname twxy ${FILENAME} xytwvaporfile.010
#tmp_num=("001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" "015" "016" "017" "018" "019" "020" "021" "022" "023" "024" "025" "026" "027" "028" "029" "030")
#tmp_num=("001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" "015" "016" "017" "018" "019" "020" "021" "022" "023" "024" "025" "026" "027" "028" "029" "030" "031" "032" "033" "034" "035" "036" "037" "038" "039" "040" "041" "042" "043" "044" "045" "046" "047" "048" "049" "050" "051" "052" "053" "054" "055" "056" "057" "058" "059" "060" "061" "062" "063" "064" "065" "066" "067" "068" "069" "070" "071" "072" "073" "074" "075" "076" "077" "078" "079" "080" "081")
#for i in {0..29}
#do
#printf ${tmp_num[$i]}\n
#raw2vdf -swapbytes -ts $i -varname bx ${FILENAME} ../../../EFLUX/d0072/B3D_R.${tmp_num[$i]}.BX
#raw2vdf -swapbytes -ts $i -varname by ${FILENAME} ../../../EFLUX/d0072/B3D_R.${tmp_num[$i]}.BY
#raw2vdf -swapbytes -ts $i -varname bz ${FILENAME} ../../../EFLUX/d0072/B3D_R.${tmp_num[$i]}.BZ
##raw2vdf -swapbytes -ts $i -varname vx ${FILENAME} ../../../EFLUX/d0072/B3D_R.${tmp_num[$i]}.VX
##raw2vdf -swapbytes -ts $i -varname vy ${FILENAME} ../../../EFLUX/d0072/B3D_R.${tmp_num[$i]}.VY
##raw2vdf -swapbytes -ts $i -varname vz ${FILENAME} ../../../EFLUX/d0072/B3D_R.${tmp_num[$i]}.VZ
#raw2vdf -swapbytes -ts $i -varname cx ${FILENAME} ../../../EFLUX/d0072/B3D_R.${tmp_num[$i]}.CX
#raw2vdf -swapbytes -ts $i -varname cy ${FILENAME} ../../../EFLUX/d0072/B3D_R.${tmp_num[$i]}.CY
#raw2vdf -swapbytes -ts $i -varname cz ${FILENAME} ../../../EFLUX/d0072/B3D_R.${tmp_num[$i]}.CZ
##raw2vdf -swapbytes -ts $i -varname pr ${FILENAME} ../../../EFLUX/d0072/B3D_R.${tmp_num[$i]}.PR
##raw2vdf -swapbytes -ts $i -varname ro ${FILENAME} ../../../EFLUX/d0072/B3D_R.${tmp_num[$i]}.RO
#raw2vdf -ts $i -varname twyz ${FILENAME} ../../tmp_python/retest_yz/3d_data_yztwist/yztw3d${tmp_num[$i]}
#raw2vdf -ts $i -varname twxy ${FILENAME} ../../tmp_python/retest_xy/xytwvaporfile.${tmp_num[$i]}
#printf "NOW TIME STEP %d\n" $i
#printf "===============\n"
#done
