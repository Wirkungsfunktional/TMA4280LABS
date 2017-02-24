CC=mpiCC


mpi_run:
	$(CC) -o main main.cpp
	mpirun -n 4 ./main

mpi_omp_run:
	$(CC) -o main_mpi_omp main.cpp -fopenmp
	mpirun -n 4 ./main_mpi_omp

