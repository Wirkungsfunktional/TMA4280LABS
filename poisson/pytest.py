#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import dst, idst









def main(args):
    n = 100
    m = n-1
    h= 1/ n
    x = np.arange(0, 1+h/2, h)
    d = 2.0 * (1.0 - np.cos((np.arange(0, m, 1) + 1) * np.pi / n))
    tt = 3
    
    b = np.zeros((m,m))
    bt = np.zeros((m,m))
    ana_u = np.zeros((m,m))
    
    for i in range(m):
        for j in range(m):
            b[i][j] = 5*np.pi**2 * np.sin(np.pi*h*i) * np.sin(2*np.pi*h*j)
            ana_u[i][j] = np.sin(np.pi*h*i) * np.sin(2*np.pi*h*j)
            
        
    for i in range(m):
        b[i] = dst(b[i], type=tt)
    bt = np.transpose(b)
    for i in range(m):
        bt[i] = idst(bt[i], type=tt)
        
        
    
    for i in range(m):
        for j in range(m):
            bt[i][j] = bt[i][j] / (d[i] + d[j]);

    for i in range(m):
        bt[i] = dst(b[i], type=tt)
    b = np.transpose(bt)
    for i in range(m):
        b[i] = idst(bt[i], type=tt)
    
    print np.sum(np.abs(b - ana_u))
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
