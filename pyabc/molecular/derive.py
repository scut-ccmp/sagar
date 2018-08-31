from pyabc.molecular.structure import Molecular
from pyabc.molecular import utils
from itertools import combinations

import numpy as np


class ConfigurationGenerator(object):
    '''
    这个类用于产生各种替换原子的需求
    '''

    def __init__(self, mol, symprec=1e-3):
        if not isinstance(mol, Molecular):
            raise TypeError(
                "want pyabc.molecular.structure.Molecular, got {:}".
                format(type(mol)))
        self._pmolecular = mol
        self._pres = symprec

    def get_configurations(self, sites, e_num, perms):
        '''
        e_num: concentration of impurity
        sites: possibility atom type of each site
        '''
        if len(e_num) == 1:
            print('There is nothing to substitute')
            return
        arg_site = np.array([len(ii) for ii in sites])
        ind = np.where(arg_site == 1)
        if not ind[0].size:  # all sites can be substituted
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
