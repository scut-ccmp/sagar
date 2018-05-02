# coding: utf-8
# Distributed under the terms of the MIT License.

import unittest
import numpy as np

from ababe.stru.scaffold import ModifiedCell
from ababe.stru.clarifier import AtomRemoveClarifier, VerboseAtomRemoveClarifier
from ababe.stru.element import Specie

class TestAtomRemoveClarifier(unittest.TestCase):

    def setUp(self):
        self.latt = np.array([[4.898979, 0.000000, 0.000000],
                              [2.449490,  4.242641, 0.000000],
                              [1.632993, -0.000000, 4.618802]])
        self.pos = np.array([[0.208333, 0.333333, 0.375000],
                             [0.375000, 0.000000, 0.875000],
                             [0.541667, 0.666667, 0.375000],
                             [0.708333, 0.333333, 0.875000],
                             [0.875000, 0.000000, 0.375000],
                             [0.000000, 0.000000, 0.000000],
                             [0.166667, 0.666667, 0.500000],
                             [0.333333, 0.333333, 0.000000],
                             [0.500000, 0.000000, 0.500000],
                             [0.666667, 0.666667, 0.000000],
                             [0.833333, 0.333333, 0.500000],
                             [0.041667, 0.666667, 0.875000]])
        self.numbers = np.array([16,16,16,16,16,30,30,30,30,30,30,55])
        self.modcell = ModifiedCell(self.latt, self.pos, self.numbers)

        centers = np.array([[0.041667, 0.666667, 0.875000]])
        ele = Specie('Zn')
        r = 2
        self.nearZnClarifier = AtomRemoveClarifier(centers, r, Specie('Zn'))

    def test_clarify(self):
        expect_pos = np.array([[0.208333, 0.333333, 0.375000],
                               [0.375000, 0.000000, 0.875000],
                               [0.541667, 0.666667, 0.375000],
                               [0.708333, 0.333333, 0.875000],
                               [0.875000, 0.000000, 0.375000],
                               [0.500000, 0.000000, 0.500000],
                               [0.833333, 0.333333, 0.500000],
                               [0.041667, 0.666667, 0.875000]])
        expect_numbers = np.array([16,16,16,16,16,30,30,55])
        expect_newcell = ModifiedCell(self.latt, expect_pos, expect_numbers)

        newcell = self.nearZnClarifier.clarify(self.modcell)
        self.assertEqual(newcell, expect_newcell)


class TestVerboseAtomRemoveClarifier(unittest.TestCase):

    def setUp(self):
        self.latt = np.array([[4.898979, 0.000000, 0.000000],
                              [2.449490,  4.242641, 0.000000],
                              [1.632993, -0.000000, 4.618802]])
        self.pos = np.array([[0.208333, 0.333333, 0.375000],
                             [0.375000, 0.000000, 0.875000],
                             [0.541667, 0.666667, 0.375000],
                             [0.708333, 0.333333, 0.875000],
                             [0.875000, 0.000000, 0.375000],
                             [0.000000, 0.000000, 0.000000],
                             [0.166667, 0.666667, 0.500000],
                             [0.333333, 0.333333, 0.000000],
                             [0.500000, 0.000000, 0.500000],
                             [0.666667, 0.666667, 0.000000],
                             [0.833333, 0.333333, 0.500000],
                             [0.041667, 0.666667, 0.875000]])
        self.numbers = np.array([16,16,16,16,55,30,30,30,30,30,30,55])
        self.modcell = ModifiedCell(self.latt, self.pos, self.numbers)

        ele = Specie('Zn')
        r = 2
        self.nearZnClarifier = VerboseAtomRemoveClarifier(Specie('Cs'), r, Specie('Zn'))

    def test_clarify(self):
        expect_pos = np.array([[0.208333, 0.333333, 0.375000],
                               [0.375000, 0.000000, 0.875000],
                               [0.541667, 0.666667, 0.375000],
                               [0.708333, 0.333333, 0.875000],
                               [0.875000, 0.000000, 0.375000],
                               [0.041667, 0.666667, 0.875000]])
        expect_numbers = np.array([16,16,16,16,55,55])
        expect_newcell = ModifiedCell(self.latt, expect_pos, expect_numbers)

        newcell = self.nearZnClarifier.clarify(self.modcell)
        self.assertEqual(newcell, expect_newcell)


if __name__ == "__main__":
    import nose2
    nose2.main()
