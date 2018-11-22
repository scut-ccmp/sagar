# -*- coding: utf-8 -*-
import unittest
import numpy

from sagar.io.vasp import _read_string, _write_string
from sagar.crystal.structure import Cell


class TestCellIO(unittest.TestCase):

    def setUp(self):
        self.string = '''H1 He1
 1.000000
    0.000000   0.500000   0.533333
    0.500000   0.000000   0.566667
    0.500000   0.500000   0.000000
H He
1 1
Direct
  0.000000 0.000000 0.000000 H
  0.250000 0.250000 0.250000 He
'''
        latt = [0., 0.5, 0.53333333333,
                0.5, 0., 0.56666666667,
                0.5, 0.5, 0.]
        positions = [0., 0., 0.,
                     0.25, 0.25, 0.25]
        atoms = ['H', 'He']
        self.cell = Cell(latt, positions, atoms)

    def test_read_string(self):
        data = self.string
        got_cell = _read_string(data)

        wanted_latt = self.cell.lattice
        wanted_positions = self.cell.positions
        wanted_atoms = self.cell.atoms

        self.assertTrue(numpy.allclose(got_cell.lattice, wanted_latt))
        self.assertTrue(numpy.allclose(got_cell.positions, wanted_positions))
        self.assertTrue(numpy.allclose(got_cell.atoms, wanted_atoms))

    def test_write_string(self):
        # wanted = self.string
        # got = _write_string(self.cell, long_format=False)
        #
        # self.assertEqual(got, wanted)
        pass

    def test_write_string_vacc(self):
        wanted = '''H1
 1.000000
    0.000000   0.500000   0.533333
    0.500000   0.000000   0.566667
    0.500000   0.500000   0.000000
H
1
Direct
  0.000000 0.000000 0.000000 H
'''
        latt = [0., 0.5, 0.53333333333,
                0.5, 0., 0.56666666667,
                0.5, 0.5, 0.]
        positions = [0., 0., 0.,
                     0.25, 0.25, 0.25]
        atoms = ['H', 'Vacc']
        cell = Cell(latt, positions, atoms)

        got = _write_string(cell, long_format=False)

        self.assertEqual(got, wanted)

    def test_Cart_string(self):
        data = '''C2
1.0
        4.9208798409         0.0000000000         0.0000000000
       -2.4604404879         4.2616066235         0.0000000000
        0.0000000000         0.0000000000        10.0000000000
    C
    8
Cartesian
     0.000000000         0.000000000         0.000000000
    -1.230220199         2.130803347         0.000000000
     2.460439920         0.000000000         0.000000000
     1.230219722         2.130803347         0.000000000
     1.230232477         0.710260868         0.000000000
     0.000012159         2.841064692         0.000000000
     3.690671921         0.710260868         0.000000000
     2.460451603         2.841064692         0.000000000
'''
        latt = [[4.9208798409 ,        0.0000000000 ,        0.0000000000],
       [-2.4604404879   ,      4.2616066235     ,    0.0000000000],
       [ 0.0000000000    ,     0.0000000000   ,     10.0000000000]]
        positions = [[0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
       [1.32703105e-08, 5.00000008e-01, 0.00000000e+00],
       [5.00000000e-01, 0.00000000e+00, 0.00000000e+00],
       [5.00000013e-01, 5.00000008e-01, 0.00000000e+00],
       [3.33335086e-01, 1.66665047e-01, 0.00000000e+00],
       [3.33335131e-01, 6.66665167e-01, 0.00000000e+00],
       [8.33334989e-01, 1.66665047e-01, 0.00000000e+00],
       [8.33335035e-01, 6.66665167e-01, 0.00000000e+00]]
        atoms = ['C' for i in range(8)]
        wanted = Cell(latt, positions, atoms)
        got = _read_string(data)
        self.assertLess(numpy.sum(abs(got.positions-wanted.positions)),1e-4)


if __name__ == "__main__":
    import nose2
    nose2.main()
