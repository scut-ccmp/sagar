import unittest
import numpy

from pyabc.io.vasp import _read_cell_from_string


class TestCellIO(unittest.TestCase):

    def test(self):
        pass

    def test_read_cell_from_string(self):
        data = """S6Zn2
1
  2.828427   0.000000   0.000000
  0.000000   2.828427   0.000000
  0.000000   0.000000  12.000000
S Zn
6 2
direct
  0.000000   0.500000   0.250000 S
  0.000000   0.500000   0.583333 S
  0.000000   0.500000   0.916667 S
  0.500000   0.000000   0.083333 S
  0.500000   0.000000   0.416667 S
  0.500000   0.000000   0.750000 S
  0.000000   0.000000   0.666667 Zn
  0.500000   0.500000   0.833333 Zn
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
        # wanted_atoms = [16, 16, 16, 16, 16, 16, 30, 30]

        self.assertTrue(numpy.allclose(got_cell.lattice, wanted_latt))
        self.assertTrue(numpy.allclose(got_cell.positions, wanted_positions))
        # self.assertTrue(np.allclose(numbers, expect_numbers))


if __name__ == "__main__":
    import nose2
    nose2.main()
