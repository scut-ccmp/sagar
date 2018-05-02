# coding: utf-8
# Distributed under the terms of the MIT License.
import unittest
import os
import numpy as np

from ababe.io.yaml import YamlInput, YamlOutput
from ababe.stru.scaffold import GeneralCell

testdata_dir = os.path.join(os.path.dirname(__file__), "test_files")

class TestYamlInput(unittest.TestCase):

    def setUp(self):
        self.filename = os.path.join(testdata_dir, "zns.yaml")

    def test_init(self):
        yaml_in = YamlInput(self.filename)

    def test_from_string(self):
        yaml_in = YamlInput(self.filename)
        content = """comment: zb_ZnS
lattice:
  - [0.0, 0.5, 0.5]
  - [0.5,  0.0,  0.5]
  - [ 0.5, 0.5,  0.0]
positions:
  - [0., 0., 0.]
  - [0.25666666, 0.25, 0.25]
numbers: [30, 16]
zoom: 4
"""
        latt, pos, numbers = yaml_in.from_string(content)
        expect_latt = np.array([[0., 2., 2.],
                                [2., 0., 2.],
                                [2., 2., 0.]])
        expect_positions = np.array([[0., 0., 0.],
                                     [0.256667, 0.25, 0.25]])
        expect_numbers = np.array([30, 16])

        self.assertTrue(np.allclose(latt, expect_latt))
        self.assertTrue(np.allclose(pos, expect_positions))
        self.assertTrue(np.allclose(numbers, expect_numbers))



    def test_get_gcell(self):
        yaml_in = YamlInput(self.filename)
        gcell = yaml_in.get_cell()

        expect_latt = np.array([[0., 2., 2.],
                                [2., 0., 2.],
                                [2., 2., 0.]])
        expect_positions = np.array([[0., 0., 0.],
                                     [0.25, 0.25, 0.25]])
        expect_numbers = np.array([30, 16])
        self.assertTrue(np.allclose(gcell.lattice, expect_latt))
        self.assertTrue(np.allclose(gcell.positions, expect_positions))
        self.assertTrue(np.allclose(gcell.numbers, expect_numbers))


class TestYamlOutput(unittest.TestCase):

    def setUp(self):
        latt = np.array([[0., 0.5, 0.53333333333],
                         [0.5, 0., 0.56666666667],
                         [0.5, 0.5, 0.]])
        positions = np.array([[0., 0., 0.],
                              [0.25, 0.25, 0.25]])
        numbers = np.array([30, 16])
        cell = GeneralCell(latt, positions, numbers)
        self.yaml_out = YamlOutput(cell)

    def test_get_string(self):
        expected_str = '''comment: S1Zn1
lattice:
- [0.0, 0.5, 0.533333]
- [0.5, 0.0, 0.566667]
- [0.5, 0.5, 0.0]
numbers: [30, 16]
positions:
- [0.0, 0.0, 0.0]
- [0.25, 0.25, 0.25]
zoom: 1
'''
        self.assertEqual(str(self.yaml_out), expected_str)

    def test_write_file(self):
        tmp_file = os.path.join(testdata_dir, "testing.yaml")
        self.yaml_out.write_file(tmp_file)

        with open(tmp_file, 'r') as testing_file:
            data = testing_file.read()

        self.assertEqual(data, str(self.yaml_out))
        os.remove(tmp_file)


if __name__ == "__main__":
    import nose2
    nose2.main()
