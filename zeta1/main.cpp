


#include <iostream>
#include <sstream>

#include "../serial_sum.hpp"
#include "../mpi_sum.hpp"



int main(int argc, char* argv[])
{

    double t1, t2;
    t1 = MPI_Wtime();

    MPI_execute<zeta> mpi_op = MPI_execute<zeta>(100000);
    mpi_op.run();
    //mpi_op.run_reduce_sum();
    t2 = MPI_Wtime();
    std::cout << t2 - t1  << "\n";

    return 0;
}
