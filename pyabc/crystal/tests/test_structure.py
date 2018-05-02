import unittest

from pyabc.crystal.structure import Cell

class testCell(unittest.TestCase):

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
        atoms = [6]*12
        self.cell = Cell(lattice, positions, atoms)

if __name__ == "__main__":
    import nose2
    nose2.main()
