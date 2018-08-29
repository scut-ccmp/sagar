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


def _check_structure(symlist, atom, seq):
    
    if len(seq) == 1:# only substitute one atom
        if seq[0] != np.size(atom):
            raise ValueError('seq[0] should be equal to %d'%(np.size(atom),))
        temp = np.sort(symlist[:, atom], axis=1)
        ind = np.where(temp[:,0] == min(temp[:,0]))
        temp = sortrows(temp[ind])
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
    
    
def gene_list(m, n):
    elelist = np.arange(0, m)
    elelist = np.reshape(elelist, (m, 1))
    for k in range(0, n-1):
        midlist1 = np.kron(np.ones((m, 1)), elelist)
        midlist2 = np.kron(np.ones((1, m)), elelist[:, 0])
        midlist2 = np.sort(midlist2).transpose()
        elelist = np.hstack((midlist2, midlist1))
    return elelist

    
    
    
    
    
    
    
