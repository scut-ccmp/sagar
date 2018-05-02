# coding: utf-8
# Distributed under the terms of the MIT License.

import unittest

from ababe.stru.element import Specie, GhostSpecie

class SpecieTestCase(unittest.TestCase):

    def setUp(self):
        self.b = Specie("B")

    def test_properties(self):
        self.assertEqual(self.b.symbol, "B")
        self.assertEqual(self.b.name, "B")
        self.assertEqual(self.b.Z, 5)
        self.assertEqual(self.b.atom_mass, 10.811)
        self.assertEqual(self.b.atom_radius, 1.17)

    def test_to_sp(self):
        self.assertEqual(Specie.from_num(5), Specie("B"))
        self.assertEqual(Specie.from_num(0), GhostSpecie())

    def test_to_name(self):
        self.assertEqual(Specie.to_name(5), "B")
        self.assertEqual(Specie.to_name(0), "G")

    def test_equal(self):
        self.assertIsNot(Specie("B"), Specie("B"))
        self.assertEqual(Specie("B"), Specie("B"))

class GhostSpecieTestCase(unittest.TestCase):

    def setUp(self):
        self.g = GhostSpecie()

    def test_properties(self):
        self.assertIsInstance(self.g, Specie)
        self.assertEqual(self.g.symbol, "G")
        self.assertEqual(self.g.Z, 0)
        self.assertEqual(self.g.atom_mass, 0)
        self.assertEqual(self.g.atom_radius, 0)

    def test_equal(self):
        self.assertEqual(GhostSpecie(), GhostSpecie())

if __name__ == "__main__":
    import nose2
    nose2.main()
