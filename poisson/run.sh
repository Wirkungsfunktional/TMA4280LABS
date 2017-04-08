#!/bin/bash

make mpi_compile
mpirun -n $1 ./main_mpi $2
hexdump -v -e '1/8 "%f" "\n" ' test.txt  > final.txt
sed -r 's/,+/./g' final.txt > final2.txt
python plot_solution.py final2.txt $2

cat final2.txt | grep -n \*
wc -l final2.txt
rm test.txt final*
