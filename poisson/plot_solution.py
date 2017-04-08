#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def main(args):
    data = np.loadtxt(args[1])
    plt.imshow(data.reshape((int(args[2])-1,int(args[2])-1)) )
    plt.show()
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
