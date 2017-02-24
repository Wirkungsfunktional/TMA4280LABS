
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

    my_sum<zeta, double> s = my_sum<zeta, double>();
    std::cout << s.eval(x) << "\n";

    return 0;
}
