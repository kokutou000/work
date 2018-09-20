#!/bin/sh
#PBS -q S
#PBS -b 16
#PBS -l elapstim_req=02:00:00
#PBS -m bea
#PBS -M n-ishiguro@isee.nagoya-u.ac.jp

#PBS -v MPIPROGINF=ALL_DETAIL

cd /S/data00/G6173/e0774/inj_eflux/20180904/sh1_exec/
mpirun -nnp 4 ../nl3d_mpi_16.256.513 > debug.dat