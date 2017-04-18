#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def main(args):
    
    data = np.loadtxt(args[1])
    N = 2**np.arange(7, 14.5, 1)
    plt.loglog(N ,data)
    plt.show()
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
