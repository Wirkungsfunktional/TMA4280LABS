#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def main(args):
    
    t1 = np.loadtxt("P1/time.txt")
    t2 = np.loadtxt("P2/time.txt")
    t4 = np.loadtxt("P4/time.txt")
    t8 = np.loadtxt("P8/time.txt")
    t12 = np.loadtxt("P12/time.txt")
    t16 = np.loadtxt("P16/time.txt")
    t20 = np.loadtxt("P20/time.txt")
    t20l = np.array( [0] + t20.tolist()  )
    t24 = np.loadtxt("P24/time.txt")
    t24l = np.array( [0] + t24.tolist()  )
    t28 = np.loadtxt("P28/time.txt")
    t28l = np.array( [0,0,0] + t28.tolist()  )
    t32 = np.loadtxt("P32/time.txt")
    t36 = np.loadtxt("P36/time.txt")
    t36l = np.array( [0,0,0] + t36.tolist()  )
    
    N = 2**np.arange(7, 14.5, 1)
    PT = np.array([t1, t2, t4, t8, t12, t16, t20l, t24l, t28l, t32, t36l])
    P = np.array([1, 2, 4, 8, 12, 16, 20, 24, 28, 32, 36])
    """

    plt.loglog(N ,t1, linewidth=4.0,label=r'$P=1$')
    plt.loglog(N ,t2, linewidth=4.0,label=r'$P=2$')
    plt.loglog(N ,t4, linewidth=4.0,label=r'$P=4$')
    plt.loglog(N ,t8, linewidth=4.0,label=r'$P=8$')
    #plt.loglog(N ,t12, linewidth=4.0,label=r'$P=1$')
    plt.loglog(N ,t16, linewidth=4.0,label=r'$P=16$')
    #plt.loglog(N[1:] ,t20, linewidth=4.0,label=r'$P=1$')
    #plt.loglog(N[1:] ,t24, linewidth=4.0,label=r'$P=1$')
    #plt.loglog(N[3:] ,t28, linewidth=4.0,label=r'$P=1$')
    plt.loglog(N ,t32, linewidth=4.0,label=r'$P=32$')
    #plt.loglog(N[3:] ,t36, linewidth=4.0,label=r'$P=1$')
    plt.loglog(N ,N**2/10**5, 'k--', linewidth=4.0,label=r'$N^2$')
    
    plt.xlabel(r'$N$', fontsize=30)
    plt.ylabel(r'$t/s$', fontsize=30)
    plt.tick_params(labelsize=25)
    plt.legend(loc='best', fontsize=30)
    plt.show()
    
    
    
    symb = [('ko', 'k--'), ('ro', 'r--'), ('bo', 'b--'), ('go', 'g--')]
    fig, ax = plt.subplots()
    for i in range(4, len(N)):
        ax.plot(P, PT[0,i]/PT[:,i], symb[i-4][0], ms=15, label=r'$N=$' + str(N[i]))
        ax.plot(P, PT[0,i]/PT[:,i], symb[i-4][1])
    
    
    
    #ax.set_xscale('log', basex=2)
    plt.xlabel(r'$P$', fontsize=30)
    plt.ylabel(r'$\frac{T_1}{T_P}$', fontsize=30)
    plt.tick_params(labelsize=25)
    plt.legend(loc='best', fontsize=30)
    plt.show()
    
    
    
    
    for i in range(4, len(N)):
        plt.plot(P, PT[0,i]/PT[:,i] / P, symb[i-4][0], ms=15, label=r'$N=$' + str(N[i]))
        plt.plot(P, PT[0,i]/PT[:,i] / P, symb[i-4][1])
        

    plt.ylim(0, 1)
    plt.xlabel(r'$P$', fontsize=30)
    plt.ylabel(r'$\frac{T_1}{P T_P}$', fontsize=30)
    plt.tick_params(labelsize=25)
    plt.legend(loc='best', fontsize=30)  
    plt.show()
    
    """
    

    P = np.array([4, 9, 18, 36])
    t = np.array([320.2261798, 164.0446110, 89.8513949, 39.5570113])
    t2 = np.array([304.8037648, 158.4780198, 87.9514769, 39.5570113])
    plt.plot(P, t, 'o', ms=15, label=r'T = 36 / P')
    #plt.plot(P, t2, 'o', ms=15, label=r'T = 1')
        

    plt.xlabel(r'$P$', fontsize=30)
    plt.ylabel(r'$t/s$', fontsize=30)
    plt.tick_params(labelsize=25)
    plt.legend(loc='best', fontsize=30)  
    plt.show()
    
    
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
