
compile:
	gcc -c poisson.c -g
	gfortran poisson.o fst.f -o main

mpi_compile:
	mpicc -c mpi_poisson.c -g
	mpif90 mpi_poisson.o fst.f -o main_mpi
