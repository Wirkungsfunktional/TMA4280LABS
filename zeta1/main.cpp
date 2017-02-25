


#include <iostream>
#include <sstream>

#include "../serial_sum.hpp"
#include "../mpi_sum.hpp"



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
    MPI_execute<zeta> mpi_op = MPI_execute<zeta>(x);
    mpi_op.run();


    return 0;
}
