#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 11:24:07 2018

@author: hecc
"""


import unittest
import numpy as np
from pyabc.molecular import utils
from pyabc.molecular import structure


class TestUtils(unittest.TestCase):

    def test_check_struc(self):
        positions = np.loadtxt('pos.txt')
        atoms = ['C'] * 60
        molecular = structure.Molecular(positions, atoms)
        all_perms = molecular.get_symmetry_permutation(0.05)
        from itertools import combinations
        all_comb = combinations(range(60), 2)
        all_type = []
        for comb in all_comb:
            all_type.append(utils._check_structure(all_perms, comb, (2,)))
        all_type = np.unique(np.array(all_type), axis=0)
        self.assertEqual(np.shape(all_type)[0], 23)


if __name__ == "__main__":
    import nose2
    nose2.main()
