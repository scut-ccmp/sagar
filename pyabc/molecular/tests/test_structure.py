#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import unittest
import numpy
from pyabc.molecular.structure import Molecular


class TestCell(unittest.TestCase):

    def setUp(self):
        positions = numpy.loadtxt('z12.txt')
        atoms = ['C'] * 20
        self.molecular = Molecular(positions, atoms)

    def test_init(self):

        positions = [0, 0, 0]
        atoms = ['NaN_10']
        mole = Molecular(positions, atoms)
        self.assertEqual(mole.atoms.tolist(), [1010])

    def test_get_permutations(self):
        pos = numpy.loadtxt('pos.txt')
        atoms = ['C'] * 60
        c60 = Molecular(pos, atoms)
        self.assertEqual(numpy.shape(c60.get_symmetry_permutation(pres=0.05))
                         [0], 120)

    def test_check(self):
        si_pos = [-0.01, -0, 0,
                  0, 0, 0]
        si_atoms = [14, 14]
        mole = Molecular(si_pos, si_atoms)
        self.assertTrue(mole.check(limit=0.1))

        si_pos = [-0.125, -0.125, -0.125,
                  0.125, 0.125, 0.125,
                  0, 0.125, 0.125]
        si_atoms = [14, 14, 8]
        mole = Molecular(si_pos, si_atoms)
        self.assertFalse(mole.check(limit=0.1))


if __name__ == "__main__":
    import nose2
    nose2.main()
