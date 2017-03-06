


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


    MPI_execute<mach, 1000> mpi_op(x, world_size, world_rank);

    double t1, t2;
    t1 = MPI_Wtime();
    mpi_op.run();
    t2 = MPI_Wtime();
    if (world_rank == 0)
        std::cout  << "Time: " << (t2 - t1) << "\n";

    MPI_Finalize();

    return 0;
}
