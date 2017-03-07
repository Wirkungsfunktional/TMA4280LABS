


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
    MPI_Init(NULL, NULL);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    MPI_execute<mach, 0> mpi_op(x, world_size, world_rank);
    mpi_op.run();

    MPI_Finalize();

    return 0;
}
