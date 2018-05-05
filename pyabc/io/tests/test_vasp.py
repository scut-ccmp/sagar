import unittest
import numpy

from pyabc.io.vasp import _read_cell_from_string


class TestCellIO(unittest.TestCase):

    def test(self):
        pass

    def test_read_cell_from_string(self):
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
        got_cell = _read_cell_from_string(data)

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


if __name__ == "__main__":
    import nose2
    nose2.main()
