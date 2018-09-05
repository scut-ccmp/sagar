#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from collections import namedtuple
from itertools import combinations


def sortrows(x):
    # TODO to be modified
    (m, n) = np.shape(x)
    ndx = np.arange(m)
    for ii in range(n-1, -1, -1):
        tmp = x[ndx, ii]
        q = np.argsort(tmp)
        ndx = ndx[q]
    x = x[ndx]
    return x

def check_structure(symlist, atom, seq):
    if len(seq) == 1:# only substitute one atom
        if seq[0] != np.size(atom):
            raise ValueError('seq[0] should be equal to %d'%(np.size(atom),))
        temp = np.sort(symlist[:, atom], axis=1)
        ind = np.where(temp[:,0] == min(temp[:,0]))
        temp = np.unique(temp[ind],axis=0)
        return temp[0]
    else:# two or more kinds of atoms
        if sum(seq) != np.size(atom):
            raise ValueError('sum of seq should be equal to %d'%(np.size(atom),))
        temp_seq = [(0, seq[0])]
        for ii in range(1, len(seq)-1):
            temp_seq.append((sum(seq[:ii])-1, sum(seq[:ii+1])))
        temp_seq.append((sum(seq[:len(seq)])-1, sum(seq)))
        all_poslist = symlist[:, atom]
        for ii in temp_seq:
            temp = range(ii[0], ii[1])
            all_poslist[:,temp] = np.sort(all_poslist[:,temp], axis=1)
        ind = np.where(all_poslist[:,0] == min(all_poslist[:,0]))
        all_poslist = sortrows(all_poslist[ind])
        return all_poslist[0]

def get_config(pos, sites, perms, e_num):
    all_type, tmp_site, e_num = get_all_type(e_num, sites, perms)
    unique_type = np.unique(np.array(all_type), axis=0)
    impurity_atom = []
    for ii in range(1, len(e_num)):
        impurity_atom.extend([tmp_site[ii]]*e_num[ii])
    all_structure = []
    tmp_mole = namedtuple('Molecule', ['pos','atoms', 'deg'])
    for single_type in unique_type:
        tmp_pos = pos[single_type]
        tmp_ind = np.setdiff1d(np.arange(len(sites)),
                               single_type)
        tmp_pos = np.vstack((tmp_pos, pos[tmp_ind]))
        tmp_atom = [tmp_site[0]]*sum(e_num)
        ind = 0
        for ii in single_type:
            tmp_atom[ii] = impurity_atom[ind]
            ind += 1
        all_structure.append((tmp_pos, tmp_atom))
    unique_type_dege = []
    indd = 0
    for single_type in unique_type:
        ind = np.where((all_type == single_type).all(axis=1))
        unique_type_dege.append(tmp_mole(all_structure[indd][0],
                                         all_structure[indd][1],
                                         len(ind[0])))
        indd += 1
    return unique_type_dege

def get_all_type(e_num, sites, perms):
    arg_site = np.array([len(ii) for ii in sites])
    e_num = np.array(e_num)
    ind = np.where(arg_site == 1)
    if not ind[0].size: # all sites can be substituted
        _ind = np.argsort(-e_num)
        e_num_new, site = [e_num[i] for i in _ind], sites[0]
        all_comb = combinations(range(sum(e_num_new)), sum(e_num_new[1:]))
        tmp_site = tuple(site[i] for i in _ind)
    else:
        need_ind = np.setdiff1d(np.arange(len(sites)), ind)
        _ind = np.argsort(-e_num)
        e_num_new, site = [e_num[i] for i in _ind], sites[need_ind[0]]
        all_comb = combinations(need_ind, sum(e_num_new[1:]))
        tmp_site = tuple(site[i] for i in _ind)
    all_type = []
    for comb in all_comb:
        all_type.append(check_structure(perms, comb,
                                               [i for i in e_num_new[1:]]))
    return all_type, tmp_site, e_num_new
