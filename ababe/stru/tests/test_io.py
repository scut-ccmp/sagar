# encoding: utf-8
# Distributed under the terms of the MIT License.

import os
import unittest

import numpy as np
from math import sqrt
from ababe.stru.scaffold import CStru
from ababe.stru.io import VaspPOSCAR, YamlOutput
from ababe.stru.scaffold import GeneralCell


class testYamlOutput(unittest.TestCase):

    def setUp(self):
        a_latt = np.array([[0.5, 0.5, -0.5],
                           [-0.5, 0.5, 0.5],
                           [0.5, -0.5, 0.5]])
        a_pos = np.array([[0, 0, 0]])
        a_num = np.array([16])
        a_cell = GeneralCell(a_latt, a_pos, a_num)
        self.a_yaml = YamlOutput(a_cell, zoom=3)

    def test_get_string(self):

        expected_str = '''comment: S1
lattice:
- [0.5, 0.5, -0.5]
- [-0.5, 0.5, 0.5]
- [0.5, -0.5, 0.5]
numbers: [16]
positions:
- [0, 0, 0]
zoom: 3
'''
        self.assertEqual(str(self.a_yaml), expected_str)

    def test_write_file(self):
        """.
        The function is tested by save structure
        to a POSCAR file, and then read from it.
        Compare the parameters read from to the
        origin input parameter. Using almostEqual
        """
        tmp_file = "YAML.testing"
        self.a_yaml.write(tmp_file)

        with open(tmp_file, 'r') as testing_file:
            data = testing_file.read()

        self.assertEqual(data, str(self.a_yaml))
        os.remove(tmp_file)


class testVaspPOSCAR(unittest.TestCase):

    def test_get_string(self):
        boron_arr = np.array([[[5, 0, 5, 5, 5, 5],
                               [5, 5, 5, 5, 5, 5],
                               [5, 5, 0, 5, 5, 5]]])
        latt = [[0, 0, 20],
                [1, 0, 0],
                [0.5, sqrt(3)/2, 0]]
        boron_stru = CStru.from_array(latt, boron_arr)
        poscar = VaspPOSCAR(boron_stru.get_gcell(), zoom=4)

        expected_str = '''B16
4
  0.000000   0.000000  20.000000
  3.000000   0.000000   0.000000
  3.000000   5.196152   0.000000
B
16
direct
  0.000000   0.000000   0.000000 B
  0.000000   0.000000   0.333333 B
  0.000000   0.000000   0.500000 B
  0.000000   0.000000   0.666667 B
  0.000000   0.000000   0.833333 B
  0.000000   0.333333   0.000000 B
  0.000000   0.333333   0.166667 B
  0.000000   0.333333   0.333333 B
  0.000000   0.333333   0.500000 B
  0.000000   0.333333   0.666667 B
  0.000000   0.333333   0.833333 B
  0.000000   0.666667   0.000000 B
  0.000000   0.666667   0.166667 B
  0.000000   0.666667   0.500000 B
  0.000000   0.666667   0.666667 B
  0.000000   0.666667   0.833333 B
'''
        self.assertEqual(str(poscar), expected_str)

        bcu_arr = np.array([[[5, 29, 5, 5, 5, 5],
                             [5, 5, 5, 5, 5, 5],
                             [5, 5, 29, 5, 5, 5]]])
        latt = [[0, 0, 20],
                [1, 0, 0],
                [0.5, sqrt(3)/2, 0]]
        bcu_stru = CStru.from_array(latt, bcu_arr)
        poscar_bcu = VaspPOSCAR(bcu_stru.get_gcell(), zoom=4)

        expected_str_bcu = '''B16Cu2
4
  0.000000   0.000000  20.000000
  3.000000   0.000000   0.000000
  3.000000   5.196152   0.000000
B Cu
16 2
direct
  0.000000   0.000000   0.000000 B
  0.000000   0.000000   0.333333 B
  0.000000   0.000000   0.500000 B
  0.000000   0.000000   0.666667 B
  0.000000   0.000000   0.833333 B
  0.000000   0.333333   0.000000 B
  0.000000   0.333333   0.166667 B
  0.000000   0.333333   0.333333 B
  0.000000   0.333333   0.500000 B
  0.000000   0.333333   0.666667 B
  0.000000   0.333333   0.833333 B
  0.000000   0.666667   0.000000 B
  0.000000   0.666667   0.166667 B
  0.000000   0.666667   0.500000 B
  0.000000   0.666667   0.666667 B
  0.000000   0.666667   0.833333 B
  0.000000   0.000000   0.166667 Cu
  0.000000   0.666667   0.333333 Cu
'''
        self.assertEqual(str(poscar_bcu), expected_str_bcu)

    def test_write_file(self):
        """.
        The function is tested by save structure
        to a POSCAR file, and then read from it.
        Compare the parameters read from to the
        origin input parameter. Using almostEqual
        """
        bcu_arr = np.array([[[5, 29, 5, 5, 5, 5],
                             [5, 5, 5, 5, 5, 5],
                             [5, 5, 29, 5, 5, 5]]])
        latt = [[0, 0, 20],
                [1, 0, 0],
                [0.5, sqrt(3)/2, 0]]
        bcu_stru = CStru.from_array(latt, bcu_arr)
        poscar_bcu = VaspPOSCAR(bcu_stru.get_gcell(), zoom=4)

        tmp_file = "POSCAR.testing"
        poscar_bcu.write(tmp_file)

        with open(tmp_file, 'r') as testing_file:
            data = testing_file.read()

        self.assertEqual(data, str(poscar_bcu))
        os.remove(tmp_file)

class TestGeneralCellInput(unittest.TestCase):

    def setUp(self):
        pass

    def test_init(self):
        pass

    def test_get_gcell(self):
        pass

if __name__ == "__main__":
    import nose2
    nose2.main()
