#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from pyabc.molecular.structure import Molecular
from pyabc.molecular import utils
from itertools import combinations

import numpy as np


class SturctureGenerator(object):
    '''
    这个类用于产生各种替换原子的需求
    '''

    def __init__(self, pmolecular, pres=1e-3):
        if not isinstance(pmolecular, Molecular):
            raise TypeError(
                "want pyabc.molecular.structure.Molecular, got {:}".
                format(type(pmolecular)))
        self._pmolecular = pmolecular
        self._pres = pres

    def _remove_redudant(self, e_num, sites, perms):
        '''
        e_num: concentration of impurity
        sites: possibility atom type of each site
        '''
        if len(e_num) == 1:
            print('There is nothing to substitute')
            return
        arg_site = np.array([len(ii) for ii in sites])
        ind = np.where(arg_site == 1)
        if not ind[0].size: # all sites can be substituted
            all_comb = combinations(range(len(sites)), sum(e_num[1:]))
        else:
            need_ind = np.setdiff1d(np.arange(len(sites)), ind)
            all_comb = combinations(need_ind, sum(e_num[1:]))
        all_type = []
        for comb in all_comb:
            all_type.append(utils._check_structure(perms, comb,
                                                   [i for i in e_num[1:]]))
        all_type = np.unique(np.array(all_type), axis=0)
        return all_type


if __name__ == "__main__":
    positions = np.loadtxt('test/pos.txt')
    atoms = ['C'] * 60
    molecular = Molecular(positions, atoms)
    perms = molecular.get_symmetry_permutation(0.05)
    S = SturctureGenerator(molecular)
    e_num = [58, 2]
    sites_1 = list([(5, 6) for i in range(20)])
    sites_2 = list([(6, ) for i in range(40)])
    sites_1.extend(sites_2)
    all_type = S._remove_redudant(e_num, sites_1, perms)
