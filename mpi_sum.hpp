#ifndef MPI_SUM_H_
#define MPI_SUM_H_


#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <cassert>






template<class FUNC, int LOOP>
class MPI_execute {

private:
    int sys_size;
    int reduced_size;
    int world_size, world_rank;

public:
    explicit MPI_execute(int N, int size, int rank){
        sys_size = N;
        world_size = size;
        world_rank = rank;
        reduced_size = sys_size / world_size;
        assert( ((world_size != 0) && ((world_size & (~world_size + 1)) == world_size)) != 0 );
    }
    void run();
    void run_reduce_sum();

};


template<class FUNC, int LOOP>
void MPI_execute<FUNC, LOOP>::run() {

    double v[reduced_size];
    double sum=0;


    if (world_rank==0) {
        double t1, t2;
        t1 = MPI_Wtime();

        for (int i=1; i<world_size; i++) {


            #pragma omp parallel for
            for (int k=0;k<reduced_size;k++) {
                v[k] = FUNC::eval( (double) (reduced_size*i + (k+1)) );
            }

            MPI_Send(v, reduced_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
        #pragma omp parallel for
        for (int k=1;k<=reduced_size;k++) {
            v[k-1] = FUNC::eval( (double) ( k ) );
        }

        for (int i=0;i<=LOOP;i++) {
            #pragma omp parallel for reduction (+:sum)
            for (int k=0;k<reduced_size;k++) {
                sum += v[k];
            }
        }
        double erg;
        for (int i=1; i<world_size; i++) {
            MPI_Recv(&erg, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum += erg;
        }
        double pi = FUNC::finalize(sum);
        t2 = MPI_Wtime();
        if (!LOOP) {
            std::cout   << "Pi: " << std::setprecision(15) << pi  << "\n"
                        << "Err: " << std::abs(M_PI - pi) << "\n"
                        << "Time: " << t2 - t1 << "\n";
        }



    } else {

        MPI_Recv(v, reduced_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i=0;i<=LOOP;i++) {
            #pragma omp parallel for reduction (+:sum)
            for (int k=0;k<reduced_size;k++) {
                sum += v[k];
            }
        }
        MPI_Send(&sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

}



template<class FUNC, int LOOP>
void MPI_execute<FUNC, LOOP>::run_reduce_sum() {

    double erg;
    double final_erg;
    double t1 = MPI_Wtime();

    for (int i=0; i<reduced_size; i++) {
        erg += FUNC::eval( (double) (reduced_size*world_rank + (i+1)) );
    }
    MPI_Reduce(&erg, &final_erg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (world_rank==0) {
        double pi = FUNC::finalize(final_erg);
        double t2 = MPI_Wtime();
        std::cout   << "Pi: " << std::setprecision(15) << pi  << "\n"
                    << "Err: " << std::abs(M_PI - pi) << "\n"
                    << "Time: " << t2 - t1 << "\n";

    }


}





#endif // MPI_SUM_H_
