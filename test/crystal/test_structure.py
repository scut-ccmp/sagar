# -*- coding: utf-8 -*-
import unittest
import numpy

from sagar.crystal.structure import Cell, MutableCell


class TestUtils(unittest.TestCase):

    def test_car_to_frac(self):
        pass

    def test_frac_to_car(self):
        pass


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
        # TODO: more complex test case needed

    def test_extend_bug(self):
        sc_latt = [4, 0, 0,
                   0, 4, 0,
                   0, 0, 4]
        sc_pos = [(0, 0, 0)]
        sc_atoms = [1]
        sc_pcell = Cell(sc_latt, sc_pos, sc_atoms)
        op_ext = numpy.array([1, 0, 0,
                              0, 1, 0,
                              0, 0, 5]).reshape((3, 3))
        ext_sc = sc_pcell.extend(op_ext)    # raise ValueError
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

    def test_check(self):
        lattice = numpy.array([0.0, 2.75, 2.75,
                               2.75, 0.0, 2.75,
                               2.75, 2.75, 0.0])
        si_pos = [-0.125, -0.125, -0.125,
                  0.125, 0.125, 0.125]
        si_atoms = [14, 14]
        cell = Cell(lattice, si_pos, si_atoms)
        self.assertTrue(cell.check(limit=0.1))

        si_pos = [-0.125, -0.125, -0.125,
                  0.125, 0.125, 0.125,
                  0.13, 0.125, 0.125]
        si_atoms = [14, 14, 8]
        cell = Cell(lattice, si_pos, si_atoms)
        self.assertFalse(cell.check(limit=0.1))
        # 测试元素
        self.assertTrue(cell.check(limit=0.1, elements=['O']))

        lattice = numpy.array([1.0, 0.0, 0.0,
                               0.0, 1.0, 0.0,
                               0.0, 0.0, 1.0])
        si_pos = [-0.01, 0.0, 0.0,
                  0.01, 0.0, 0.0]
        si_atoms = [14, 14]
        # 两点之间的距离为0.02<0.05， 但在没有考虑边界条件时距离是0.98
        cell = Cell(lattice, si_pos, si_atoms)
        self.assertFalse(cell.check(limit=0.05))


class TestMutableCell(unittest.TestCase):

    def setUp(self):
        # 使用金刚石Si作为测试模型
        self.lattice = numpy.array([0.0, 0.5, 0.5,
                                    0.5, 0.0, 0.5,
                                    0.5, 0.5, 0.0]).reshape((3, 3))

    def siteEqual(self, site_a, site_b, prec=1e-5):
        p_a, e_a = numpy.array(site_a[0]), site_a[1]
        p_b, e_b = numpy.array(site_b[0]), site_b[1]
        if e_a == e_b and numpy.linalg.norm(p_a-p_b) < prec:
            is_eq = True
        else:
            is_eq = False
        return is_eq

    def test_init(self):
        # 创建空胞
        lattice = numpy.copy(self.lattice)
        mcell = MutableCell(lattice)
        self.assertEqual(len(mcell._sites), 0)

        # 创建金刚石Si原胞
        lattice = numpy.copy(self.lattice)
        si_sites = [[(-0.125, -0.125, -0.125), "Si"],
                    [(0.125, 0.125, 0.125), "Si"]]
        mcell = MutableCell(lattice, sites=si_sites)
        self.assertEqual(len(mcell._sites), 2)

    def test_to_cell(self):
        lattice = numpy.copy(self.lattice)
        si_sites = [[(-0.125, -0.125, -0.125), "Si"],
                    [(0.125, 0.125, 0.125), "Si"]]
        mcell = MutableCell(lattice, sites=si_sites)
        c = mcell.to_cell()
        self.assertTrue(isinstance(c, Cell))
        numpy.testing.assert_almost_equal(c.positions, numpy.array([0.875, 0.875, 0.875,
                                                                    0.125, 0.125, 0.125]).reshape((2, 3)))
        numpy.testing.assert_almost_equal(c.atoms, numpy.array([14, 14]))

    def test_from_cell(self):
        lattice = numpy.copy(self.lattice)
        positions = numpy.array([0.0, 0.0, 0.0,
                                 0.25, 0.25, 0.25]).reshape((2, 3))
        atoms = numpy.array([6, 6])
        origin_lattice = numpy.copy(lattice)
        origin_positions = numpy.copy(positions)
        origin_atoms = numpy.copy(atoms)

        cell = Cell(lattice, positions, atoms)
        mcell = MutableCell.from_cell(cell)
        self.assertTrue(isinstance(mcell, MutableCell))

        # Make sure mcell modified not change the cell
        mcell.remove_site(0)
        numpy.testing.assert_almost_equal(cell.lattice, origin_lattice)
        numpy.testing.assert_almost_equal(cell.positions, origin_positions)
        numpy.testing.assert_almost_equal(cell.atoms, origin_atoms)

    def test_add_site(self):
        lattice = numpy.copy(self.lattice)
        mcell = MutableCell(lattice)
        mcell.add_site([(-0.125, -0.125, -0.125), "Zn"])
        mcell.add_site([(0.125, 0.125, 0.125), "S"])

        self.assertEqual(mcell._sites[0], [(-0.125, -0.125, -0.125), "Zn"])
        self.assertEqual(mcell._sites[1], [(0.125, 0.125, 0.125), "S"])

    def test_set_site(self):
        lattice = numpy.copy(self.lattice)
        si_sites = [[(-0.125, -0.125, -0.125), "Si"],
                    [(0.125, 0.125, 0.125), "Si"]]
        mcell = MutableCell(lattice, sites=si_sites)
        mcell.set_site(0, [(-0.125, -0.125, -0.125), "Zn"])
        mcell.set_site(1, [(0.126, 0.125, 0.125), "S"])

        self.assertEqual(len(mcell._sites), 2)
        self.assertEqual(mcell._sites[0], [(-0.125, -0.125, -0.125), "Zn"])
        self.assertEqual(mcell._sites[1], [(0.126, 0.125, 0.125), "S"])

    def test_remove_site(self):
        lattice = numpy.copy(self.lattice)
        si_sites = [[(-0.125, -0.125, -0.125), "Si"],
                    [(0.125, 0.125, 0.125), "Si"]]
        mcell = MutableCell(lattice, sites=si_sites)
        mcell.remove_site(1)
        self.assertEqual(len(mcell._sites), 2)
        self.assertEqual(mcell._sites[0], [(-0.125, -0.125, -0.125), "Si"])
        self.assertEqual(mcell._sites[1], [(0.125, 0.125, 0.125), "Vacc"])

    def test_rotate_site_by_z(self):
        lattice = numpy.copy(self.lattice)
        si_sites = [[(-0.125, -0.125, -0.125), "Si"],
                    [(0.125, 0.125, 0.125), "Si"]]
        mcell = MutableCell(lattice, sites=si_sites)
        mcell.rotate_site_by_z(1, (0,0,0), 90)
        self.assertTrue(self.siteEqual(mcell._sites[1], [(0.125, -0.125, 0.125), "Si"]))

        mcell = MutableCell(lattice, sites=si_sites)
        mcell.rotate_site_by_z(1, (0.125,0,0), -90)
        self.assertTrue(self.siteEqual(mcell._sites[1], [(0.25, 0.0, 0.125), "Si"]))

        mcell = MutableCell(lattice, sites=[[(0.6, 0.6, 0.5), "C"]])
        mcell.rotate_site_by_z(0, (0.5, 0.5,0), -90)
        self.assertTrue(self.siteEqual(mcell._sites[0], [(0.4, 0.6, 0.5), "C"]))

    def test_get_site(self):
        lattice = numpy.copy(self.lattice)
        si_sites = [[(-0.125, -0.125, -0.125), "Si"],
                    [(0.125, 0.125, 0.125), "Si"]]
        mcell = MutableCell(lattice, sites=si_sites)
        self.assertEqual(mcell.get_site(0), [(-0.125, -0.125, -0.125), "Si"])
        self.assertEqual(mcell.get_site(1), [(0.125, 0.125, 0.125), "Si"])

    def test_get_car_site(self):
        lattice = numpy.copy(numpy.array([0.0, 1.5, 1.5,
                                          1.5, 0.0, 1.5,
                                          1.5, 1.5, 0.0]).reshape(3, 3))
        si_sites = [[(-0.125, -0.125, -0.125), "Si"],
                    [(0.125, 0.125, 0.125), "Si"]]
        mcell = MutableCell(lattice, sites=si_sites)
        self.assertEqual(mcell.get_car_site(0), [(-0.375, -0.375, -0.375), "Si"])
        self.assertEqual(mcell.get_car_site(1), [(0.375, 0.375, 0.375), "Si"])


if __name__ == "__main__":
    import nose2
    nose2.main()
