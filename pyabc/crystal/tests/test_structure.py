import unittest

from pyabc.crystal.structure import Cell
# from pyabc.crystal.structure import is_primitive_cell, get_primitive_cell


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

    def test_is_primitive(self):
        bcc_latt = [0.5, 0.5, -0.5,
                    -0.5, 0.5, 0.5,
                    0.5, -0.5, 0.5]
        bcc_pos = [(0, 0, 0)]
        bcc_atoms = [0]
        bcc_pcell = Cell(bcc_latt, bcc_pos, bcc_atoms)
        self.assertTrue(bcc_pcell.is_primitive())

    def test_get_primitive(self):
        fcc_latt = [5, 0, 0,
                    0, 5, 0,
                    0, 0, 5]
        fcc_pos = [(0, 0, 0),
                   (0, 0.5, 0.5),
                   (0.5, 0, 0.5),
                   (0.5, 0.5, 0)]
        fcc_atoms = [0, 0, 0, 0]
        con_cell = Cell(fcc_latt, fcc_pos, fcc_atoms)
        pcell = con_cell.get_primitive_cell()
        self.assertEqual(pcell.atoms, [0])


if __name__ == "__main__":
    import nose2
    nose2.main()
