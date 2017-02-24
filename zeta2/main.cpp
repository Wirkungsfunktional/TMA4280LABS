


#include <iostream>
#include <sstream>

#include "../serial_sum.hpp"
#include "../omp_sum.hpp"



int main(int argc, char* argv[])
{
    OMP_execute<zeta> omp_op = OMP_execute<zeta>(1000);
    omp_op.run();

    return 0;
}
