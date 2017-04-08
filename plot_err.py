#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
from matplotlib import pyplot as plt

err_zeta = np.array([
0.40298,
0.21898,
0.114295,
0.0583996,
0.0295188,
0.0148399,
0.00744013,
0.00372513,
0.00186383,
0.000932232,
0.000466195,
0.000233117,
0.000116564,
5.8283e-05,
2.91418e-05,
1.4571e-05,
7.28552e-06,
3.64276e-06,
1.82138e-06,
9.10692e-07,
4.55346e-07,
2.27673e-07,
1.13837e-07,
5.69228e-08])

N = np.array([
2,
4,
8,
16,
32,
64,
128,
256,
512,
1024,
2048,
4096,
8192,
16384,
32768,
65536,
131072,
262144,
524288,
1048576,
2097152,
4194304,
8388608,
16777216])

err_mach = np.array([
0.000995624,
8.81408e-07,
1.19016e-12,
1.33227e-15,
1.33227e-15,
1.33227e-15,
1.33227e-15,
1.33227e-15,
1.33227e-15,
1.33227e-15,
1.33227e-15,
1.33227e-15,
1.33227e-15,
1.33227e-15,
1.33227e-15,
1.33227e-15,
1.33227e-15,
1.33227e-15,
1.33227e-15,
1.33227e-15,
1.33227e-15,
1.33227e-15,
1.33227e-15,
1.33227e-15
])






def main(args):
    plt.loglog(N, err_zeta, linewidth=4.0, label='zeta')
    plt.loglog(N, err_mach, linewidth=4.0, label='mach')
    
    plt.ylabel(r'$| \pi_n - \pi |$', fontsize=30)
    plt.xlabel(r'$n$', fontsize=30)
    plt.tick_params(labelsize=25)
    plt.legend(loc='best', fontsize=30)
    plt.show()
    
    
    
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
