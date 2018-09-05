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
        # check input parameters
        if len(self._pmolecular.atoms) != sum(e_num):
            raise ValueError('Atom number of molecular %s are not equal\
                             to e_num %s'%(len(self._pmolecular.atoms), 
                             sum(e_num)))
        if len(e_num) == 1:
            print('There is nothing to substitute')
            return
        arg_site = np.array([len(ii) for ii in sites])
        if np.unique(arg_site).size > 2:
            raise ValueError('the `len` of tuples in sites can not be \
                             three or more')
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
    
    def get_degeneracy(self, e_num, sites, perms):
        '''
        to get degeneracy of the configuration in all combinations
        
        e_num: concentration of impurity
        sites: possibility atom type of each site
        
        return list including tuples, for each tuple:
        tuple[0] is the configuration
        tuple[1] is the degeneracy
        '''
                # check input parameters
        if len(self._pmolecular.atoms) != sum(e_num):
            raise ValueError('Atom number of molecular %s are not equal\
                             to e_num %s'%(len(self._pmolecular.atoms), 
                             sum(e_num)))
        if len(e_num) == 1:
            print('There is nothing to substitute')
            return
        arg_site = np.array([len(ii) for ii in sites])
        if np.unique(arg_site).size > 2:
            raise ValueError('the `len` of tuples in sites can not be \
                             three or more')
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
        unique_type = np.unique(np.array(all_type), axis=0)
        unique_type_dege = []
        for single_type in unique_type:
            ind = np.where((all_type == single_type).all(axis=1))
            unique_type_dege.append((single_type, len(ind[0])))
        return unique_type_dege
        
        
if __name__ == "__main__":
    positions = np.loadtxt('test/pos.txt')
    atoms = ['C'] * 60
    molecular = Molecular(positions, atoms)
    perms = molecular.get_symmetry_permutation(0.05)
    S = SturctureGenerator(molecular)
    e_num = [57, 3]
    sites_1 = list([(5, 6) for i in range(60)])
    type_dege = S.get_degeneracy(e_num, sites_1, perms)
