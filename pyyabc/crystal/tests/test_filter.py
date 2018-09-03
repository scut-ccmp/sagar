import unittest
import numpy

from pyyabc.crystal.filter import *
from pyyabc.crystal.structure import Cell
from pyyabc.crystal.mutate import cell_to_mcell


class TestMinDistanceFilter(unittest.TestCase):

    def setUp(self):
        lattice = numpy.array([7.1, 0, 0,
                               0, 7.1, 0,
                               0, 0, 7.1]).reshape((3, 3))
        positions = numpy.array([0.00, 0.00, 0.00,
                                 0.00, 0.00, 0.50,
                                 0.00, 0.50, 0.00,
                                 0.00, 0.50, 0.50,
                                 0.50, 0.00, 0.00,
                                 0.50, 0.00, 0.50,
                                 0.50, 0.50, 0.00,
                                 0.50, 0.50, 0.50,
                                 0.00, 0.25, 0.25,
                                 0.00, 0.25, 0.75,
                                 0.00, 0.75, 0.25,
                                 0.00, 0.75, 0.75,
                                 0.50, 0.25, 0.25,
                                 0.50, 0.25, 0.75,
                                 0.50, 0.75, 0.25,
                                 0.50, 0.75, 0.75,
                                 0.25, 0.00, 0.25,
                                 0.25, 0.00, 0.75,
                                 0.25, 0.50, 0.25,
                                 0.25, 0.50, 0.75,
                                 0.75, 0.00, 0.25,
                                 0.75, 0.00, 0.75,
                                 0.75, 0.50, 0.25,
                                 0.75, 0.50, 0.75,
                                 0.25, 0.25, 0.00,
                                 0.25, 0.25, 0.50,
                                 0.25, 0.75, 0.00,
                                 0.25, 0.75, 0.50,
                                 0.75, 0.25, 0.00,
                                 0.75, 0.25, 0.50,
                                 0.75, 0.75, 0.00,
                                 0.75, 0.75, 0.50])
        atoms = numpy.array([28] * 32)
        cell = Cell(lattice, positions, atoms)

        mcell = cell_to_mcell(cell)
        mcell.set_site(0, ([0, 0, 0], 'Co'))
        mcell.set_site(16, ([0.25, 0, 0.25], 'Co'))
        self.cell_nearest_1 = mcell.to_cell()

        mcell = cell_to_mcell(cell)
        mcell.set_site(0, ([0, 0, 0], 'Co'))
        mcell.set_site(21, ([0.75, 0, 0.75], 'Co'))
        self.cell_nearest_2 = mcell.to_cell()

        mcell = cell_to_mcell(cell)
        mcell.set_site(1, ([0, 0, 0.5], 'Co'))
        mcell.set_site(5, ([0.5, 0, 0.5], 'Co'))
        self.cell_far = mcell.to_cell()

    def test_is_accepted(self):
        filter = MinDistanceFilter('Co', 2.6)

        # 胞内和原点最近
        self.assertFalse(filter.is_accepted(self.cell_nearest_1))

        # 胞外和原点最近
        self.assertFalse(filter.is_accepted(self.cell_nearest_2))

        # 不相邻,分散
        self.assertTrue(filter.is_accepted(self.cell_far))


class TestSpaceGroupFilter(unittest.TestCase):

    def setUp(self):
        lattice = numpy.array([5.5, 0, 0,
                               0, 5.5, 0,
                               0, 0, 5.5]).reshape((3, 3))
        positions = numpy.array([0.0, 0.0, 0.0,
                                 0.0, 0.5, 0.5,
                                 0.5, 0.0, 0.5,
                                 0.5, 0.5, 0.0,
                                 0.5, 0.5, 0.5,
                                 0.5, 0.0, 0.0,
                                 0.0, 0.5, 0.0,
                                 0.0, 0.0, 0.5])
        atoms = numpy.array([11, 11, 11, 11, 17, 17, 17, 17])
        self.cell = Cell(lattice, positions, atoms)

    def test_is_accepted(self):
        true_spg_filter = SpaceGroupFilter([225, ], symprec=1e-3)
        self.assertTrue(true_spg_filter.is_accepted(self.cell))

        false_spg_filter = SpaceGroupFilter([0, 1, 227, ], symprec=1e-3)
        self.assertFalse(false_spg_filter.is_accepted(self.cell))
