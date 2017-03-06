CC=mpiCC


mpi_run:
	$(CC) -o main main.cpp -O2


mpi_omp_run:
	$(CC) -o main_mpi_omp main.cpp -fopenmp

mpi_loop:
	$(CC) -o main_loop main_loop.cpp -O2
	for i in 1 2 4; do \
      echo $$i;\
      mpirun -n $$i ./main_loop 1000000;\
  done
