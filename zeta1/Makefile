CC=mpiCC


mpi_run:
	$(CC) -o main main.cpp


mpi_omp_run:
	$(CC) -o main_mpi_omp main.cpp -fopenmp

mpi_loop:
	for i in 1 2 4 8 16 32; do \
      echo $$i;\
			mpirun -n $$i ./main 10000000;\
  done
