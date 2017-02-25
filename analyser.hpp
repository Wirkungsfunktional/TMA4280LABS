#ifndef ANALYSER_H_
#define ANALYSER_H_

#include "serial_sum.hpp"
#include "mpi_sum.hpp"
#include "omp_sum.hpp"




struct omp_wrapper
{
    static inline void run(int N) {
        OMP_execute<zeta> omp_op = OMP_execute<zeta>(N);
        omp_op.run();
    }
};






template<class FUNCTOR>
class analyser {

private:
    int ID;

public:
    explicit analyser() {
        ID = 0;
    }

    void run();
};


template<class FUNCTOR>
void analyser<FUNCTOR>::run() {
    int n;
    for (int i=1;i<=24; i++) {
        n = 1 << i;
        std::cout << i << " " << n << "\n";
        FUNCTOR::run( n );
    }
}







#endif //ANALYSER_H_
