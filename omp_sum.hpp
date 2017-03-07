#ifndef OMP_SUM_H_
#define OMP_SUM_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <omp.h>
#include <iostream>
#include <iomanip>

template<class FUNC>
class OMP_execute {

private:
    int sys_size;

public:
    explicit OMP_execute(int N){
        sys_size = N;
    }
    void run();
};

template<class FUNC>
void OMP_execute<FUNC>::run() {
    double ta, te;
    ta = omp_get_wtime();
    double v[sys_size];

    #pragma omp parallel for
    for (int i=1; i<=sys_size; i++) {
        v[i-1] = FUNC::eval( (double) i);
    }

    double erg = 0.0;

    #pragma omp parallel for reduction (+:erg)
    for (int i=0; i<sys_size; i++) {
        erg += v[i];
    }
    double pi = FUNC::finalize(erg);
    te = omp_get_wtime();
    std::cout    << "Pi1: " << std::setprecision(15)<< pi  << "\n"
                << "Err: " << std::abs(M_PI - pi) << "\n"
                << "Time: " << te - ta << "\n";
}


#endif // OMP_SUM_H_
