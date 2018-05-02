# encoding: utf-8
# Distributed under the terms of the MIT License.

import os
import unittest

import numpy as np
from math import sqrt
from ababe.stru.scaffold import CStru
from ababe.stru.scaffold import GeneralCell
from ababe.io.io import GeneralIO
from ababe.io.vasp import VaspOutput
from ababe.io.yaml import YamlOutput

testdata_dir = os.path.join(os.path.dirname(__file__), "test_files")

class TestGeneralIO(unittest.TestCase):

    def setUp(self):
        latt = np.array([[0., 0.5, 0.53333333333], [0.5, 0., 0.56666666667],
                        [0.5, 0.5, 0.]])
        pos = np.array([[0., 0., 0.],
                             [0.25, 0.25, 0.25]])
        numbers = np.array([30, 16])
        gcell = GeneralCell(latt, pos, numbers)
        self.vasp_out = VaspOutput(gcell)
        self.yaml_out = YamlOutput(gcell)
        self.io_s = GeneralIO(gcell)

    def test_vasp_from_file(self):
        filename = os.path.join(testdata_dir, "POSCAR_")
        gcell = GeneralIO.from_file(filename)
        expect_latt = np.array([[2.8284, 0.0000, 0.0000],
                                [0.0000, 2.8284, 0.0000],
                                [0.0000, 0.0000, 12.0000]])
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

    def test_yaml_from_file(self):
        filename = os.path.join(testdata_dir, "zns.yaml")
        gcell = GeneralIO.from_file(filename)

        expect_latt = np.array([[0., 2., 2.],
                                [2., 0., 2.],
                                [2., 2., 0.]])
        expect_positions = np.array([[0., 0., 0.],
                                     [0.25, 0.25, 0.25]])
        expect_numbers = np.array([30, 16])
        self.assertTrue(np.allclose(gcell.lattice, expect_latt))
        self.assertTrue(np.allclose(gcell.positions, expect_positions))
        self.assertTrue(np.allclose(gcell.numbers, expect_numbers))

    def test_vasp_write_file(self):
        ##################################################
        # vasp POSCAR write test
        ##################################################
        tmp_file = os.path.join(testdata_dir, "testing.vasp")
        self.io_s.write_file(tmp_file)

        with open(tmp_file, 'r') as testing_file:
            data = testing_file.read()

        self.assertEqual(data, str(self.vasp_out))
        os.remove(tmp_file)

        #################
        tmp_file = os.path.join(testdata_dir, "POSCAR______")
        self.io_s.write_file(tmp_file)

        with open(tmp_file, 'r') as testing_file:
            data = testing_file.read()

        self.assertEqual(data, str(self.vasp_out))
        os.remove(tmp_file)

        ###################
        tmp_file = os.path.join(testdata_dir, "i_have_noname")
        self.io_s.write_file(tmp_file, fmt='vasp')

        with open(tmp_file, 'r') as testing_file:
            data = testing_file.read()

        self.assertEqual(data, str(self.vasp_out))
        os.remove(tmp_file)

    def test_yaml_write_file(self):
        ##################################################
        # yaml write test
        ##################################################
        tmp_file = os.path.join(testdata_dir, "testing.yaml")
        self.io_s.write_file(tmp_file)

        with open(tmp_file, 'r') as testing_file:
            data = testing_file.read()

        self.assertEqual(data, str(self.yaml_out))
        os.remove(tmp_file)

        ########################3
        tmp_file = os.path.join(testdata_dir, "testing")
        self.io_s.write_file(tmp_file, fmt='yaml')

        with open(tmp_file, 'r') as testing_file:
            data = testing_file.read()

        self.assertEqual(data, str(self.yaml_out))
        os.remove(tmp_file)

    def test_noname_write_file(self):
        self.io_s.write_file(fmt="vasp")
        pass


if __name__ == "__main__":
    import nose2
    nose2.main()
