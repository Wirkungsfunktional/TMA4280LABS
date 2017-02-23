
#include <iostream>
#include <sstream>

#include "../serial_sum.hpp"





int main(int argc, char** argv)
{
    std::istringstream ss(argv[1]);
    int x;
    if (!(ss >> x))
        std::cerr << "Invalid number " << argv[1] << '\n';

    sum<mach, double> s = sum<mach, double>();
    std::cout << s.eval(x) << "\n";
    s.verification_test();


    return 0;
}
