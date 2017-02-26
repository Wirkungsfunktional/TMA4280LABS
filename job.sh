#!/bin/bash


#PBS -N NAME
#PBS -A imf_tma4280
#PBS -l walltime 00:01:00
#PBS -l select=2:ncpus=32:mpiprocs=16

cd $PBS_O_WORKDIR
load openmpi

mpirun ./main 10000
