#!/bin/bash -f

#PBS -q NODEX1
#PBS -l elapstim_req=72:00:00
##PBS -l memsz_job=120GB
#PBS -T intmpi
##PBS -v OMP_STACKSIZE=512m
##PBS -v OMP_NUM_THREADS=2
##--- このOMPなんちゃらはOpenMPのやつなので
##--- MPIを利用している現状では関係ない、はず。。

##--- 計算実行、終了時に連絡をくれるらしいが
##--- いらなければコメントアウトでおｋ
##PBS -M n-ishiguro@isee.nagoya-u.ac.jp
##PBS -m be
##PBS -j o

ulimit -s unlimited
cd $PBS_O_WORKDIR
##mpirun ${NQSII_MPIOPTS} -n 4 ./a.out < input.txt > output
##mpirun ${NQSII_MPIOPTS} -n 4 ./nl3d_mpi02k_064x064x032
mpirun ${NQSII_MPIOPTS} -n 16 ../nl3d_mpi02k_512x512x256 > debug.dat
