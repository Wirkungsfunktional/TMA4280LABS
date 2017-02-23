
#include <iostream>

#include "../serial_sum.hpp"






int main(int argc, char* argv[])
{
    //sum<zeta, double> s = sum<zeta, double>();
    //std::cout << s.eval(1000) << "\n";
    //s.verification_test();

    MPI_execute<zeta> mpi_op = MPI_execute<zeta>(1000);
    //mpi_op.run();
    mpi_op.run_reduce_sum();

    return 0;
}
