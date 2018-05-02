# coding: utf-8
# Distributed under the terms of the MIT License.
import unittest
import os
import numpy as np

from ababe.io.vasp import VaspInput, VaspOutput
from ababe.stru.scaffold import GeneralCell

testdata_dir = os.path.join(os.path.dirname(__file__), "test_files")

class TestVaspInput(unittest.TestCase):

    def setUp(self):
        self.filename = os.path.join(testdata_dir, "POSCAR_")

    def test_init(self):
        vasp_in = VaspInput(self.filename)

    def test_from_string(self):
        vasp_in = VaspInput(self.filename)
        content = """S6Zn2
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
        latt, pos, numbers = vasp_in.from_string(content)
        expect_latt = np.array([[2.828427, 0.000000, 0.000000],
                                [0.000000, 2.828427, 0.000000],
                                [0.000000, 0.000000, 12.000000]])
        expect_positions = np.array([[0.000000,   0.500000,   0.250000],
                                     [0.000000,   0.500000,   0.583333],
                                     [0.000000,   0.500000,   0.916667],
                                     [0.500000,   0.000000,   0.083333],
                                     [0.500000,   0.000000,   0.416667],
                                     [0.500000,   0.000000,   0.750000],
                                     [0.000000,   0.000000,   0.666667],
                                     [0.500000,   0.500000,   0.833333]])
        expect_numbers = [16, 16, 16, 16, 16, 16, 30, 30]

        self.assertTrue(np.allclose(latt, expect_latt))
        self.assertTrue(np.allclose(pos, expect_positions))
        self.assertTrue(np.allclose(numbers, expect_numbers))

    def test_get_gcell(self):
        vasp_in = VaspInput(self.filename)
        gcell = vasp_in.get_cell()
        expect_latt = np.array([[2.828427, 0.000000, 0.000000],
                                [0.000000, 2.828427, 0.000000],
                                [0.000000, 0.000000, 12.000000]])
        expect_positions = np.array([[ 0.    ,  0.    ,  0.666667],
                                     [ 0.    ,  0.5   ,  0.250000],
                                     [ 0.    ,  0.5   ,  0.583333],
                                     [ 0.    ,  0.5   ,  0.916667],
                                     [ 0.5   ,  0.    ,  0.083333],
                                     [ 0.5   ,  0.    ,  0.416667],
                                     [ 0.5   ,  0.    ,  0.750000],
                                     [ 0.5   ,  0.5   ,  0.833333]])

        expect_numbers = [30, 16, 16, 16, 16, 16, 16, 30]

        self.assertTrue(np.allclose(gcell.lattice, expect_latt))
        self.assertTrue(np.allclose(gcell.positions, expect_positions))
        self.assertTrue(np.allclose(gcell.numbers, expect_numbers))


class TestVaspOutput(unittest.TestCase):

    def setUp(self):
        latt = np.array([[0., 0.5, 0.53333333333],
                         [0.5, 0., 0.56666666667],
                         [0.5, 0.5, 0.]])
        positions = np.array([[0., 0., 0.],
                              [0.25, 0.25, 0.25]])
        numbers = np.array([30, 16])
        cell = GeneralCell(latt, positions, numbers)
        self.vasp_out = VaspOutput(cell)

    def test_get_string(self):
        expected_str = '''S1Zn1
1
  0.000000   0.500000   0.533333
  0.500000   0.000000   0.566667
  0.500000   0.500000   0.000000
S Zn
1 1
direct
  0.250000   0.250000   0.250000 S
  0.000000   0.000000   0.000000 Zn
'''
        self.assertEqual(str(self.vasp_out), expected_str)

    def test_write_file(self):
        tmp_file = os.path.join(testdata_dir, "testing.vasp")
        self.vasp_out.write_file(tmp_file)

        with open(tmp_file, 'r') as testing_file:
            data = testing_file.read()

        self.assertEqual(data, str(self.vasp_out))
        os.remove(tmp_file)

if __name__ == "__main__":
    import nose2
    nose2.main()
