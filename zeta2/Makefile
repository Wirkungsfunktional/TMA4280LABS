CC=mpiCC

omp_run:
	$(CC) -o main main.cpp -fopenmp
	./main

omp_loop: 
	$(CC) -o main main.cpp -fopenmp
	for i in 1 2 4; do \
      		echo $$i;\
		export OMP_NUM_THREADS=$$i;\
		./main 1000000;\
  	done


