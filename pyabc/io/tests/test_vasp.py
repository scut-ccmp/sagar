import unittest
import numpy

from pyabc.io.vasp import _read_string, _write_string
from pyabc.crystal.structure import Cell


class TestCellIO(unittest.TestCase):

    def test(self):
        pass

    def test_read_string(self):
        data = """H6He2
1
  2.828427   0.000000   0.000000
  0.000000   2.828427   0.000000
  0.000000   0.000000  12.000000
H He
6 2
direct
  0.000000   0.500000   0.250000 H
  0.000000   0.500000   0.583333 H
  0.000000   0.500000   0.916667 H
  0.500000   0.000000   0.083333 H
  0.500000   0.000000   0.416667 H
  0.500000   0.000000   0.750000 H
  0.000000   0.000000   0.666667 He
  0.500000   0.500000   0.833333 He
"""
        got_cell = _read_string(data)

        wanted_latt = numpy.array([[2.828427, 0.000000, 0.000000],
                                   [0.000000, 2.828427, 0.000000],
                                   [0.000000, 0.000000, 12.000000]])
        wanted_positions = numpy.array([[0.000000,   0.500000,   0.250000],
                                        [0.000000,   0.500000,   0.583333],
                                        [0.000000,   0.500000,   0.916667],
                                        [0.500000,   0.000000,   0.083333],
                                        [0.500000,   0.000000,   0.416667],
                                        [0.500000,   0.000000,   0.750000],
                                        [0.000000,   0.000000,   0.666667],
                                        [0.500000,   0.500000,   0.833333]])
        wanted_atoms = numpy.array([1, 1, 1, 1, 1, 1, 2, 2])

        self.assertTrue(numpy.allclose(got_cell.lattice, wanted_latt))
        self.assertTrue(numpy.allclose(got_cell.positions, wanted_positions))
        self.assertTrue(numpy.allclose(got_cell.atoms, wanted_atoms))

    def test_write_string(self):
        wanted = '''H1 He1
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
        cell = Cell(latt, positions, atoms)
        got = _write_string(cell, long_format=False)

        self.assertEqual(got, wanted)


if __name__ == "__main__":
    import nose2
    nose2.main()
