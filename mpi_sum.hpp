#ifndef MPI_SUM_H_
#define MPI_SUM_H_


#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <mpi.h>
#include <omp.h>


template<class FUNC>
class MPI_execute {

private:
    int sys_size;

public:
    explicit MPI_execute(int N){
        sys_size = N;
    }
    void run();
    void run_reduce_sum();
    void run_reduce_sum_rec();

};


template<class FUNC>
void MPI_execute<FUNC>::run() {

    MPI_Init(NULL, NULL);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int reduced_size = sys_size / world_size;

    if (world_rank==0) {

        for (int i=1; i<world_size; i++) {
            double v[reduced_size];
            for (int k=0;k<reduced_size;k++) {
                v[k] = FUNC::eval( (double) (reduced_size*i + (k+1)) );
            }
            MPI_Send(v, reduced_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
        double sum = 0;

        for (int k=1;k<=reduced_size;k++) {
            sum += FUNC::eval( (double) ( k ) );
        }
        for (int i=1; i<world_size; i++) {
            double erg;
            MPI_Recv(&erg, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum += erg;
        }
        std::cout << FUNC::finalize(sum) << "\n";



    } else {
        double v[reduced_size];
        double sum=0;
        MPI_Recv(v, reduced_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int k=0;k<reduced_size;k++) {
            sum += v[k];
        }
        MPI_Send(&sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Finalize();
}

double rec(int N, int rank, double z, int pivot) {
    if (N==1) {
        return z;
    }
    if (rank < pivot) {
        double x = rec(N/2, rank, z, pivot/2);
        double erg;
        std::cout << rank << " send to " << rank + N/2 << "\n";
        MPI_Send(&x, 1, MPI_DOUBLE, rank + N/2, 0, MPI_COMM_WORLD);
        MPI_Recv(&erg, 1, MPI_DOUBLE, rank + N/2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        return x + erg;
    } else {
        double x = rec(N/2, rank, z, pivot + pivot/2 - 1);
        double erg;
        std::cout << rank << " send to " << rank - N/2 << "\n";
        MPI_Send(&x, 1, MPI_DOUBLE, rank - N/2, 0, MPI_COMM_WORLD);
        MPI_Recv(&erg, 1, MPI_DOUBLE, rank - N/2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        return x + erg;
    }
}
/*
double collect(int N, double x) {
    int deepth = get_depth(N);
    double erg=x;
    double z;
    for (int i=1; i<=deepth; i++) {
        if (rank % 2 == 0) {
            MPI_Send(&erg, 1, MPI_DOUBLE, rank + pow(2, i-1), 0, MPI_COMM_WORLD);
            MPI_Recv(&z, 1, MPI_DOUBLE, rank + pow(2, i-1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            erg += z;
        } else {
            MPI_Send(&erg, 1, MPI_DOUBLE, rank - pow(2, i-1), 0, MPI_COMM_WORLD);
            MPI_Recv(&z, 1, MPI_DOUBLE, rank - pow(2, i-1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            erg += z;
        }
    }

}*/




template<class FUNC>
void MPI_execute<FUNC>::run_reduce_sum() {

    MPI::Init();
    int size = MPI::COMM_WORLD.Get_size();
    int rank = MPI::COMM_WORLD.Get_rank();
    int reduced_size = sys_size / size;
    double erg;
    double final_erg;

    for (int i=0; i<reduced_size; i++) {
        erg += FUNC::eval( (double) (reduced_size*rank + (i+1)) );
    }
    //MPI::COMM_WORLD.Reduce(&erg, &final_erg, 1, MPI::DOUBLE, MPI::SUM, 0);
    final_erg = rec(size, rank, erg, size/2);
    //if (rank==0) {
    std::cout << "Rank " << rank << ": " << FUNC::finalize(final_erg) << "\n";
    //}
    MPI::Finalize();
}





#endif // MPI_SUM_H_
