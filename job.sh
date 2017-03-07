#!/bin/bash


#PBS -N NAME
#PBS -A imf_tma4280
#PBS -l walltime=00:01:00
#PBS -l select=1:ncpus=32:mpiprocs=16

cd $PBS_O_WORKDIR

module load gcc
module load openmpi

mpirun ./main 10000
