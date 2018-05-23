import unittest
import numpy

from pyabc.crystal.structure import Cell


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
        fcc_latt = [0, 5, 5,
                    5, 0, 5,
                    5, 5, 0]
        fcc_pos = [(0, 0, 0)]
        fcc_atoms = [0]
        fcc_pcell = Cell(fcc_latt, fcc_pos, fcc_atoms)
        # op_ext = numpy.array([-1, 1, 1,
        #                       1, -1, 1,
        #                       1, 1, -1]).reshape((3, 3))
        op_ext = numpy.array([1, 0, 1,
                              0, 2, 1,
                              0, 0, 2]).reshape((3, 3))
        ext_fcc = fcc_pcell.extend(op_ext)

        wanted_latt = numpy.array([5, 10, 5,
                                   15, 5, 10,
                                   10, 10, 0]).reshape((3, 3))
        wanted_pos = numpy.array([(0, 0, 0),
                                  (0, 0, 0.5),
                                  (0, 0.5, 0.25),
                                  (0, 0.5, 0.75)])
        wanted_atoms = numpy.array([0, 0, 0, 0])
        numpy.testing.assert_almost_equal(ext_fcc.lattice, wanted_latt)
        numpy.testing.assert_almost_equal(ext_fcc.positions, wanted_pos)
        numpy.testing.assert_almost_equal(ext_fcc.atoms, wanted_atoms)

    def test_extend_not_hnf(self):
        fcc_latt = [0, 5, 5,
                    5, 0, 5,
                    5, 5, 0]
        fcc_pos = [(0, 0, 0),
                   (0.5, 0.5, 0.5)]
        fcc_atoms = [1, 2]
        fcc_pcell = Cell(fcc_latt, fcc_pos, fcc_atoms)
        op_ext = numpy.array([-1, 1, 1,
                              1, -1, 1,
                              1, 1, -1]).reshape((3, 3))
        ext_fcc = fcc_pcell.extend(op_ext)

        wanted_latt = numpy.array([10, 0, 0,
                                   0, 10, 0,
                                   0, 0, 10]).reshape((3, 3))
        wanted_pos = numpy.array([(0, 0, 0),
                                  (0.5, 0.5, 0),
                                  (0.5, 0, 0.5),
                                  (0, 0.5, 0.5),
                                  (0.5, 0.5, 0.5),
                                  (0, 0, 0.5),
                                  (0, 0.5, 0),
                                  (0.5, 0, 0)])
        wanted_atoms = numpy.array([1, 1, 1, 1, 2, 2, 2, 2])
        numpy.testing.assert_almost_equal(ext_fcc.lattice, wanted_latt)
        numpy.testing.assert_almost_equal(ext_fcc.positions, wanted_pos)
        numpy.testing.assert_almost_equal(ext_fcc.atoms, wanted_atoms)

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
