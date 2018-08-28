#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 26 15:47:24 2018

@author: hecc
"""


import numpy as np


def gene_list(m, n):
    elelist = np.arange(0, m)
    elelist = np.reshape(elelist, (m, 1))
    for k in range(0, n-1):
        midlist1 = np.kron(np.ones((m, 1)), elelist)
        midlist2 = np.kron(np.ones((1, m)), elelist[:, 0])
        midlist2 = np.sort(midlist2).transpose()
        elelist = np.hstack((midlist2, midlist1))
    return elelist


def sortrows(a):
    arr = []
    (row, col) = np.shape(a)
    for ii in range(row):
        arr.append(tuple(a[ii]))
    key = [(str(i), int) for i in range(col)]
    arr = np.array(arr, dtype=key)
    # 这个是实现从一列最高优先, 其次是第二列， 第三列，直到最后一列
    arr = np.sort(arr, order=[str(i) for i in range(np.shape(a)[1])])
    return arr


def check_structure(symlist, atom, seq):
    temp_seq = [(0, seq[0])]
    for ii in range(1, len(seq)):
        temp_seq.append((seq[ii-1]+1, seq[ii]))
    all_poslist = symlist[:, atom]
    for ii in range(len(seq)):
        temp = range(temp_seq[ii][0], temp_seq[ii][1]+1)
        all_poslist[temp, :] = np.sort(all_poslist[temp], axis=1)
    all_poslist = sortrows(all_poslist)
    return all_poslist[0]


if __name__ == "__main__":
    import structure
    from itertools import combinations
    pos = np.loadtxt('z12.txt')
    atoms = ['C']*20
    C = structure.Cell(pos, atoms)
    symlist = C.get_symmetry_permutation()
    all_perm = combinations(range(20), 2)
    all_type = []
    for atom in all_perm:
        seq = [2]
        all_type.append(check_structure(symlist, atom, seq))
    all_type = np.unique(all_type, axis=0)
    print(np.shape(all_type)[0])
