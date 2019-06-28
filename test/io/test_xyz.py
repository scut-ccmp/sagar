# -*- coding: utf-8 -*-
import unittest
import numpy

from sagar.io.xyz import _read_string, _write_string
from sagar.molecule.structure import Molecule

class TestMolIO(unittest.TestCase):

    def setUp(self):
        self.string = '''2
H1He1
H   0.000000 0.000000 0.000000
He  0.250000 0.250000 0.250000
'''
        positions = [0., 0., 0.,
                     0.25, 0.25, 0.25]
        atoms = ['H', 'He']
        self.mol = Molecule(positions, atoms)

    def test_read_string(self):
        data = self.string
        got_mol = _read_string(data)

        wanted_positions = self.mol.positions
        wanted_atoms = self.mol.atoms

        self.assertTrue(numpy.allclose(got_mol.positions, wanted_positions))
        self.assertTrue(numpy.allclose(got_mol.atoms, wanted_atoms))

    def test_write_string(self):
        wanted = '''1
H1
 H 0.000000 0.000000 0.000000
'''
        positions = [0., 0., 0.,
                     2.5, 2.5, 2.5]
        atoms = ['H', 'Vac']
        mol = Molecule(positions, atoms)

        got = _write_string(mol, long_format=False)

        self.assertEqual(got, wanted)
