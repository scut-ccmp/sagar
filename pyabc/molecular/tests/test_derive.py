#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import unittest
import numpy
from pyabc.molecular import structure
from pyabc.molecular.derive import SturctureGenerator


class TestDerive(unittest.TestCase):

    def setUp(self):
        positions = numpy.loadtxt('pos.txt')
        atoms = ['C'] * 60
        self.molecular = structure.Molecular(positions, atoms)
        self.SturctureGenerator = SturctureGenerator(self.molecular, pres=1e-4)

    def test_remove_redundant_1(self):
        positions = numpy.loadtxt('pos.txt')
        atoms = ['C'] * 60
        molecular = structure.Molecular(positions, atoms)
        perms = molecular.get_symmetry_permutation(0.05)
        S = SturctureGenerator(molecular)
        e_num = [58,2]
        sites = [(5, 6) for i in range(60)]
        all_type = S._remove_redudant(e_num, sites, perms)
        self.assertEqual(numpy.shape(all_type)[0], 23)

    def test_remove_redundant_2(self):
        positions = numpy.loadtxt('pos.txt')
        atoms = ['C'] * 60
        molecular = structure.Molecular(positions, atoms)
        perms = molecular.get_symmetry_permutation(0.05)
        S = SturctureGenerator(molecular)
        e_num = [58, 2, 1]
        sites = [(5, 6, 7) for i in range(60)]
        all_type = S._remove_redudant(e_num, sites, perms)
        self.assertEqual(numpy.shape(all_type)[0], 871)

    def test_remove_redundant_3(self):
        positions = numpy.loadtxt('pos.txt')
        atoms = ['C'] * 60
        molecular = structure.Molecular(positions, atoms)
        perms = molecular.get_symmetry_permutation(0.05)
        S = SturctureGenerator(molecular)
        e_num = [58, 2]
        sites_1 = list([(5, 6) for i in range(20)])
        sites_2 = list([(6, ) for i in range(40)])
        sites_1.extend(sites_2)
        all_type = S._remove_redudant(e_num, sites_1, perms)
        self.assertEqual(numpy.shape(all_type)[0], 17)


if __name__ == "__main__":
    import nose2
    nose2.main()
