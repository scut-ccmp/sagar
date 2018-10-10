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



if __name__ == "__main__":
    import nose2
    nose2.main()
