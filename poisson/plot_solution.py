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
5.404161e-08
])

N = np.array([
64, 
128, 
256, 
512, 
1024 
])
h = 1/N

def main(args):
    plt.loglog(h, err)
    plt.plot(h, h*h)
    plt.show()
    
    #data = np.loadtxt(args[1])
    #plt.imshow(data.reshape((int(args[2])-1,int(args[2])-1)) )
    #plt.show()
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
