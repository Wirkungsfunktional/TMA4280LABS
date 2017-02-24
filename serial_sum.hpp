#ifndef SERIAL_SUM_H_
#define SERIAL_SUM_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <omp.h>


struct zeta
{
    static inline double eval(double x) {
        return 1/(x*x);
    }
    static inline double finalize(double x) {
        return sqrt(6.0*x);
    }
};

struct mach
{
    static inline double eval(double x) {
        double i = 2*x-1;
        return std::pow(-1.0, x-1) / i * (4/std::pow(5.0, i) - 1/std::pow(239.0, i));
    }
    static inline double finalize(double x) {
        return 4.0*x;
    }
};








template <class FUNC, typename T>
class my_sum {

public:
    explicit my_sum() {
        erg = T();
    }
    T eval(int n);
    T unit_test();
    void verification_test();

private:
    T erg;
};


template <class FUNC, typename T>
T my_sum<FUNC, T>::eval(int n) {
    erg = T();
    for (int i = 1; i<=n; i++) {
        erg += FUNC::eval( (T) i);
    }
    return FUNC::finalize(erg);
}

template <class FUNC, typename T>
T my_sum<FUNC, T>::unit_test() {
    return eval(3);
}

template <class FUNC, typename T>
void my_sum<FUNC, T>::verification_test() {
    T err[24];
    for (int i=1; i<=24; i++) {
        err[i] = M_PI - eval(1 << i);
        std::cout << err[i] << "\n";
    }
}













#endif // SERIAL_SUM_H_
