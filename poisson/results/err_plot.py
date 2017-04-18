#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def main(args):
    
    t1 = np.loadtxt("P1/err.txt")
    t2 = np.loadtxt("P2/err.txt")
    t4 = np.loadtxt("P4/err.txt")
    t8 = np.loadtxt("P8/err.txt")
    t12 = np.loadtxt("P12/err.txt")
    t16 = np.loadtxt("P16/err.txt")
    t20 = np.loadtxt("P20/err.txt")
    t20l = np.array( [0] + t20.tolist()  )
    t24 = np.loadtxt("P24/err.txt")
    t24l = np.array( [0] + t24.tolist()  )
    t28 = np.loadtxt("P28/err.txt")
    t28l = np.array( [0,0,0] + t28.tolist()  )
    t32 = np.loadtxt("P32/err.txt")
    t36 = np.loadtxt("P36/err.txt")
    t36l = np.array( [0,0,0] + t36.tolist()  )
    
    N = 2**np.arange(7, 14.5, 1)
    h = 1/(N+1)
    

    plt.loglog(h ,t1, linewidth=4.0,label=r'$P=1$')
    plt.loglog(h ,t2, linewidth=4.0,label=r'$P=2$')
    plt.loglog(h ,t4, linewidth=4.0,label=r'$P=4$')
    plt.loglog(h ,t8, linewidth=4.0,label=r'$P=8$')
    #plt.loglog(N ,t12, linewidth=4.0,label=r'$P=1$')
    plt.loglog(h ,t16, linewidth=4.0,label=r'$P=16$')
    #plt.loglog(N[1:] ,t20, linewidth=4.0,label=r'$P=1$')
    #plt.loglog(N[1:] ,t24, linewidth=4.0,label=r'$P=1$')
    #plt.loglog(N[3:] ,t28, linewidth=4.0,label=r'$P=1$')
    plt.loglog(h ,t32, linewidth=4.0,label=r'$P=32$')
    #plt.loglog(N[3:] ,t36, linewidth=4.0,label=r'$P=1$')
    plt.loglog(h ,h**2, 'k--', linewidth=4.0,label=r'$h^2$')
    
    plt.xlabel(r'$h$', fontsize=30)
    plt.ylabel(r'$err$', fontsize=30)
    plt.tick_params(labelsize=25)
    plt.legend(loc='best', fontsize=30)
    plt.show()
    
    
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
