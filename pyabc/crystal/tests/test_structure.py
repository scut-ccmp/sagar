import unittest

from pyabc.crystal.structure import Cell, is_primitive_cell


class TestCell(unittest.TestCase):

    def setUp(self):
        lattice = [3.0, 0, 0, 0, 2.0, 0, 0, 0, 1.0]
        positions = [(0.00000, 0.00000, 0.00000),
                     (0.00000, 0.50000, 0.00000),
                     (0.33333, 0.00000, 0.00000),
                     (0.33333, 0.50000, 0.00000),
                     (0.66666, 0.00000, 0.00000),
                     (0.66666, 0.50000, 0.00000),
                     (0.16666, 0.25000, 0.50000),
                     (0.16666, 0.75000, 0.50000),
                     (0.50000, 0.25000, 0.50000),
                     (0.50000, 0.75000, 0.50000),
                     (0.83333, 0.25000, 0.50000),
                     (0.83333, 0.75000, 0.50000)]
        atoms = ['C'] * 12
        self.cell = Cell(lattice, positions, atoms)

    def test_init(self):
        lattice = [1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0]
        positions = [0, 0, 0]
        atoms = ['NaN_10']
        cell = Cell(lattice, positions, atoms)
        self.assertEqual(cell.atoms.tolist(), [1010])

    def test_extend(self):
        pass

    def test_rotate(self):
        pass

    def test_trans(self):
        pass


class TestUtils(unittest.TestCase):

    def test_is_primitive_cell(self):
        bcc_latt = [0.5, 0.5, -0.5,
                    -0.5, 0.5, 0.5,
                    0.5, -0.5, 0.5]
        bcc_pos = [(0, 0, 0)]
        bcc_atoms = [0]
        bcc_pcell = Cell(bcc_latt, bcc_pos, bcc_atoms)
        self.assertTrue(is_primitive_cell(bcc_pcell))


if __name__ == "__main__":
    import nose2
    nose2.main()
