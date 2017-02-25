
#include <iostream>
#include <sstream>

#include "../serial_sum.hpp"






int main(int argc, char* argv[])
{
    my_sum<mach, double> s = my_sum<mach, double>();
    s.unit_test();

    return 0;
}
