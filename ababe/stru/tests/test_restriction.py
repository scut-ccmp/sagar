# coding: utf-8
# Distributed under the terms of the MIT License.

import unittest
from ababe.stru.restriction import MinDistanceRestriction
from ababe.stru.restriction import SpaceGroupRestriction
from ababe.stru.element import Specie
from ababe.stru.scaffold import GeneralCell
import numpy as np

class testMinDistanceRestriction(unittest.TestCase):

    def setUp(self):
        arr_lat = np.array([[3.0, 0, 0], [0, 2.0, 0.0], [0, 0, 1.0]])
        positions = [
                        [0.00000, 0.00000, 0.00000],
                        [0.00000, 0.50000, 0.00000],
                        [0.33333, 0.00000, 0.00000],
                        [0.33333, 0.50000, 0.00000],
                        [0.66666, 0.00000, 0.00000],
                        [0.66666, 0.50000, 0.00000],
                        [0.16666, 0.25000, 0.50000],
                        [0.16666, 0.75000, 0.50000],
                        [0.50000, 0.25000, 0.50000],
                        [0.50000, 0.75000, 0.50000],
                        [0.83333, 0.25000, 0.50000],
                        [0.83333, 0.75000, 0.50000]
                    ]
        arr_positions = np.array(positions)
        arr_numbers = np.array([5, 6, 6, 6, 5, 6,
                                6, 6, 6, 6, 6, 6])
        self.cell = GeneralCell(arr_lat, arr_positions, arr_numbers)

        self.carbon_restrc = MinDistanceRestriction((Specie('C'), 0.85))
        self.boron_restrc = MinDistanceRestriction((Specie('B'), 1.01))

    def test_is_satisfied(self):
        self.assertTrue(self.carbon_restrc.is_satisfied(self.cell))
        self.assertFalse(self.boron_restrc.is_satisfied(self.cell))


class testSpaceGroupRestriction(unittest.TestCase):

    def setUp(self):
        arr_lat = np.array([[3.0, 0, 0], [0, 2.0, 0.0], [0, 0, 1.0]])
        positions = [
                        [0.00000, 0.00000, 0.00000],
                        [0.00000, 0.50000, 0.00000],
                        [0.33333, 0.00000, 0.00000],
                        [0.33333, 0.50000, 0.00000],
                        [0.66666, 0.00000, 0.00000],
                        [0.66666, 0.50000, 0.00000],
                        [0.16666, 0.25000, 0.50000],
                        [0.16666, 0.75000, 0.50000],
                        [0.50000, 0.25000, 0.50000],
                        [0.50000, 0.75000, 0.50000],
                        [0.83333, 0.25000, 0.50000],
                        [0.83333, 0.75000, 0.50000]
                    ]
        arr_positions = np.array(positions)
        arr_numbers = np.array([5, 6, 22, 6, 5, 6,
                                6, 29, 6, 6, 6, 6])
        self.p6cell = GeneralCell(arr_lat, arr_positions, arr_numbers)

        arr_numbers = np.array([5, 6, 6, 6, 5, 6,
                                6, 6, 6, 6, 6, 6])
        self.p6orsymcell = GeneralCell(arr_lat, arr_positions, arr_numbers)

        self.p6_restrc = SpaceGroupRestriction([6])
        self.p6orsym_restrc = SpaceGroupRestriction([6, 47])

    def test_is_satisfied(self):
        self.assertTrue(self.p6_restrc.is_satisfied(self.p6cell))
        self.assertTrue(self.p6orsym_restrc.is_satisfied(self.p6orsymcell))


if __name__ == "__main__":
    import nose2
    nose2.main()
