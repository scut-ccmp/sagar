import unittest
import numpy

from pyabc.crystal.structure import Cell
# from pyabc.io.ase import read, write


class TestCellIO(unittest.TestCase):

    def setUp(self):
        self.latt = numpy.array([[4, 0, 0],
                                 [0, 4, 0],
                                 [0, 0, 4]])
        self.pos = numpy.array([(0, 0, 0),
                                (0, 0.5, 0.5)])
        self.numbers = [2, 3]

    # TODO: deadloop import, maybe ase bugs to be report.

    # def test_read(self):
    #     from ase import Atoms
    #
    #     ase_atoms = Atoms(cell=self.latt,
    #                       scaled_positions=self.pos,
    #                       numbers=self.numbers,
    #                       pbc=True)
    #     pyabc_cell = read(ase_atoms)
    #     self.assertTrue(numpy.allclose(pyabc_cell.lattice), self.latt)
    #     self.assertTrue(numpy.allclose(pyabc_cell.positions), self.pos)
    #     self.assertTrue(numpy.allclose(pyabc_cell.atoms),
    #                     numpy.array(self.numbers))
    #
    # def test_write(self):
    #     from pyabc.io.ase import write
    #     pyabc_cell = Cell(self.latt, self.pos, self.numbers)
    #     ase_atoms = write(pyabc_cell)
    #     self.assertTrue(numpy.allclose(ase_atoms.get_cell()), self.latt)
    #     self.assertTrue(numpy.allclose(
    #         ase_atoms.get_positions(wrap=True)), self.pos)
    #     self.assertTrue(numpy.allclose(pyabc_cell.atoms),
    #                     numpy.array(self.numbers))
