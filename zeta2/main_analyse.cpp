

#include <iostream>
#include <sstream>


#include "../analyser.hpp"



int main(int argc, char* argv[])
{
    analyser<omp_wrapper> ana = analyser<omp_wrapper>();
    ana.run();
    
    return 0;
}
