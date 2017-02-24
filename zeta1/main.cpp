


#include <iostream>
#include <sstream>

#include "../serial_sum.hpp"
#include "../mpi_sum.hpp"



int main(int argc, char* argv[])
{
    MPI_execute<zeta> mpi_op = MPI_execute<zeta>(1000);
    mpi_op.run();
    //mpi_op.run_reduce_sum();

    return 0;
}
