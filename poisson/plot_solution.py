#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt



err = np.array([
1.383983e-05,
3.458984e-06,
8.646850e-07,
2.161674e-07,
5.404161e-08,
1.351035e-08,
3.377277e-09,
8.470690e-10,
2.290498e-10
])

N = np.array([
64, 
128, 
256, 
512, 
1024,
2048,
4096,
8192,
8192*2
])
h = 1/N

def main(args):
    #plt.loglog(h, err, 'o')
    #plt.plot(h, h*h)
    #plt.show()
    
    data = np.loadtxt(args[1])
    plt.imshow(data.reshape((int(args[2])-1,int(args[2])-1)).transpose(),extent=[0,1,1,0])
    
    
    plt.xlabel(r'$x$', fontsize=30)
    plt.ylabel(r'$y$', fontsize=30)
    plt.tick_params(labelsize=25)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=30) 
    plt.show()
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
