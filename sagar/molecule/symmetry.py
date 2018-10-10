# -*- coding: utf-8 -*-
"""
Created on Sat Aug 25 16:02:21 2018

@author: hecc
"""
import numpy as np
from itertools import combinations, permutations


def get_distance_matrix(pos):
    n = np.shape(pos)[0]
    d = np.zeros((n, n))
    for ii in range(n):
        for jj in range(ii + 1, n):
            d[ii, jj] = np.linalg.norm(pos[ii] - pos[jj])
    d = d + d.T
    return d


def get_three_permutation(d, atoms, pres=1e-3):
    n = np.shape(d)[0]
    all_comb = combinations(range(n), 3)
    all_tri = []
    corre_comb = []
    for comb in all_comb:
        temp = [d[comb[0], comb[1]], d[comb[0], comb[2]], d[comb[1], comb[2]]]
        all_tri.append(sorted(temp))
        corre_comb.append(comb)
    all_tri = np.array(all_tri).reshape(-1, 3)
    corre_comb = np.array(corre_comb).reshape(-1, 3)
    n_all_tri = np.shape(all_tri)[0]
    temp_d = np.linalg.norm(
        all_tri - np.tile(all_tri[0], (n_all_tri, 1)), axis=1)
    ind = np.where(temp_d < pres)
    poss_com = corre_comb[ind]
    # here we have known that `1 2 3` atoms can transform into ind atoms
    # but we can not know precisely which atom is correspondant to the other
    all_perms = []
    for sub_com in poss_com:
        temp = list(permutations(sub_com))
        all_perms.append(temp)
    all_perms = np.array(all_perms).reshape(-1, 3)
    # based on the sub distance matrix
    sub_d_matrix = d[0:3, 0:3]
    tri_symm = []
    for perm in all_perms:
        dd = (sub_d_matrix[0, 1] - d[perm[0], perm[1]])**2 + \
             (sub_d_matrix[0, 2] - d[perm[0], perm[2]])**2 + \
             (sub_d_matrix[1, 2] - d[perm[1], perm[2]])**2
        if dd < pres and atoms[:3] == [atoms[perm[0]],
                                       atoms[perm[1]], atoms[perm[2]]]:
            tri_symm.append(perm)
    tri_symm = np.array(tri_symm).reshape(-1, 3)
    return tri_symm


def get_new_symm(d, origin_d, tri_symm, atoms, pres=1e-3):
    new_tri_symm = []
    n = np.shape(d)[0]
    for each_symm in tri_symm:
        left_index = np.setdiff1d(range(n), each_symm)
        for index in left_index:
            nn = np.shape(origin_d)[0]
            new_d = np.zeros((nn + 1, nn + 1))
            new_d[0:nn, 0:nn] = origin_d
            new_d[nn, 0:nn] = d[index, each_symm]
            new_d[0:nn, nn] = d[each_symm, index]
            temp_dis = np.linalg.norm(new_d - d[0:nn + 1, 0:nn + 1])
            if temp_dis < pres and atoms[nn] == atoms[index]:
                new_tri_symm.append(np.hstack((each_symm, index)))
                break
    return new_tri_symm


def get_permutations(pos, atoms, pres=1e-3):
    '''
    This algorithm are mainly derived by Chaobin Qiu
    '''
    d = get_distance_matrix(pos)
    tri_symm = get_three_permutation(d, atoms, pres=pres)
    n = np.shape(d)[0]
    for ii in range(n - 3):
        origin_d = d[0:ii + 3, 0:ii + 3]
        tri_symm = get_new_symm(d, origin_d, tri_symm, atoms, pres=pres)
    return np.array(tri_symm)
