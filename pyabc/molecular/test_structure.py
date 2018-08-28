#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 25 21:14:48 2018

@author: hecc
"""

import unittest
import numpy
from structure import Cell


class TestCell(unittest.TestCase):

    def setUp(self):
        positions = numpy.loadtxt('z12.txt')
        atoms = ['C'] * 20
        self.cell = Cell(positions, atoms)

    def test_init(self):

        positions = [0, 0, 0]
        atoms = ['NaN_10']
        cell = Cell(positions, atoms)
        self.assertEqual(cell.atoms.tolist(), [1010])

    def test_check(self):
        si_pos = [-0.01, -0, 0,
                  0, 0, 0]
        si_atoms = [14, 14]
        cell = Cell(si_pos, si_atoms)
        self.assertTrue(cell.check(limit=0.1))

        si_pos = [-0.125, -0.125, -0.125,
                  0.125, 0.125, 0.125,
                  0, 0.125, 0.125]
        si_atoms = [14, 14, 8]
        cell = Cell(si_pos, si_atoms)
        self.assertFalse(cell.check(limit=0.1))


if __name__ == "__main__":
    import nose2
    nose2.main()
