#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
def sortrows(x):
    (m, n) = np.shape(x)
    ndx = np.arange(m)
    for ii in range(n-1, -1, -1):
        q = np.argsort(x[ndx, ii])
        ndx = ndx[q]
    x = x[ndx]
    return x


if __name__ == "__main__":
    x = np.loadtxt('all_perms.txt',dtype=int)
    y = sortrows(x)