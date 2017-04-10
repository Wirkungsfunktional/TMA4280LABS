#!/bin/bash

make mpi_compile
for i in 64 128 256 512 1024; do
    echo "$i "
    mpirun -n 4 ./main_mpi $i 1
done
