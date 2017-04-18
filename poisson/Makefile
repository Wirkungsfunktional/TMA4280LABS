
compile:
	gcc -c poisson.c -g
	gfortran poisson.o fst.f -o main

mpi_compile:
	mpicc -c -fopenmp mpi_poisson.c -g -O2
	mpif90 -fopenmp mpi_poisson.o fst.f -o main_mpi
