
#include <iostream>
#include <sstream>

#include "../serial_sum.hpp"






int main(int argc, char* argv[])
{
    int x;
    if (argc < 2) {
        std::cout << "Insert number: " << "\n";
        std::cin >> x;
    } else {
        std::istringstream ss(argv[1]);
        if (!(ss >> x))
            std::cerr << "Invalid number " << argv[1] << '\n';
    }

    my_sum<mach, double> s = my_sum<mach, double>();
    std::cout << s.eval(x) << "\n";

    //MPI_execute<zeta> mpi_op = MPI_execute<zeta>(1000);
    //mpi_op.run();
    //mpi_op.run_reduce_sum();

    //OMP_execute<zeta> omp_op = OMP_execute<zeta>(1000);
    //omp_op.run();


    return 0;
}
