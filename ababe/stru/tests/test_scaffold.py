# coding: utf-8
# Distributed under the terms of the MIT License.

from ababe.stru.scaffold import SitesGrid, CStru, GeneralCell, \
                                ModifiedCell
from ababe.stru.site import Site
from ababe.stru.element import GhostSpecie, Specie
import numpy as np

import unittest

class testSitesGrid(unittest.TestCase):

    def setUp(self):
        g = GhostSpecie()
        b = Specie("B")
        self.sites   = [[[b, b],
                    [g, g]],

                   [[b, g],
                    [g, b]]]
        self.sg = SitesGrid(self.sites)
        self.allg0 = SitesGrid.sea(2, 2, 2, GhostSpecie())
        self.allg = SitesGrid.sea(4, 2, 2, GhostSpecie())

    def test_properties(self):
        # ok_(np.allclose(self.g.base_vector,
                        # np.array(m, dtype=np.float64).reshape((3,3))))
        self.assertEqual(self.allg.depth, 4)
        self.assertEqual(self.allg.length, 2)
        self.assertEqual(self.allg.width, 2)
        self.assertEqual(self.sg.sites, self.sites)

    def test_get(self):
        self.assertEqual(self.sg[0, 1, 1], GhostSpecie())

        self.sg[0, 1, 1] = Specie("B")
        self.assertEqual(self.sg[0, 1, 1], Specie("B"))

    def test_equal(self):
        other = SitesGrid.sea(4, 2, 2, Specie("Ti"))
        self.assertNotEqual(other, self.allg)
        self.assertNotEqual(None, self.allg)

    def test_deepCopy(self):
        g_copy = self.allg.deepCopy()
        self.assertEqual(g_copy, self.allg)
        self.assertNotEqual(id(g_copy), id(self.allg))
        # assert_not_equal(id(g_copy[0, 0, 0]), id(g_copy[0, 0, 0]))

    def test_hash(self):
        pass

    def test_to_array(self):
        a = np.array([0]*8).reshape((2,2,2))
        self.assertTrue(np.allclose(self.allg0.to_array(), a))
        b = np.array([5,5,0,0,5,0,0,5]).reshape((2,2,2))
        self.assertTrue(np.allclose(self.sg.to_array(), b))

    def test_from_array(self):
        arr = np.array([5,5,0,0,5,0,0,5]).reshape([2,2,2])
        ss = SitesGrid.from_array(arr)
        self.assertEqual(ss, self.sg)

    def test_random_fill(self):
        r = SitesGrid.random_fill(GhostSpecie(), (2,2,2), Specie("B"))
        self.assertIn(r.to_array().sum(), [5*x for x in range(9)])

    def test_gen_speckle(self):
        c = Specie("Cu")
        t = Specie("Ti")
        sites   = [[[c, c],
                    [t, t]],

                   [[c, t],
                    [t, c]]]
        self.sg = SitesGrid(sites)
        gen = SitesGrid.gen_speckle(Specie("Cu"), (2,2,2), Specie("Ti"), 4)
        from collections import Iterator
        self.assertIsInstance(gen, Iterator)
        self.assertIn(self.sg, gen)
        self.assertEqual(next(gen).to_array().sum(), 204)
        self.assertEqual(next(gen).to_array().sum(), 204)

class testCStru(unittest.TestCase):

    def setUp(self):
        self.m = [[1,0,0],[0,1,0],[0,0,1]]
        g = GhostSpecie()
        b = Specie("B")
        self.sites = [[[g, b],
                  [g, g]],

                 [[b, b],
                  [b, g]]]

        self.arr = np.array([0,5,0,0,5,5,5,0]).reshape([2,2,2])
        self.sg = SitesGrid(self.sites)
        self.s = CStru(self.m, self.sg)

    def test_get_property(self):
        self.assertEqual(self.s.m, self.m)
        # eq_(self.s.depth, 2)
        # eq_(self.s.width, 2)
        # eq_(self.s.length, 2)

        # eq_(self.s.get_grid, self.sites)
        # arr = [[[0, 5],
        #         [0, 0]],

        #        [[5, 5],
        #         [5, 0]]]
        # eq_(self.s.get_array(), arr)

    def test_equal(self):
        m_0 = [[1,1,1], [0,0,1], [1,0,0]]
        g = GhostSpecie()
        b = Specie("B")
        sites_0 = [[[b, b],
                    [g, g]],

                   [[b, g],
                    [g, b]]]
        sg_0 = SitesGrid(sites_0)

        diff_m = CStru(m_0, self.sg)
        diff_s = CStru(self.m, sg_0)
        self.assertEqual(self.s, self.s)
        self.assertNotEqual(diff_m, self.s)
        self.assertNotEqual(diff_s, self.s)

    def test_from_array(self):
        ss = CStru.from_array(self.m, self.arr)
        self.assertEqual(ss, self.s)

    def test_get_array(self):
        self.assertTrue(np.allclose(self.s.get_array(), self.arr))

    def test_gen_speckle(self):
        c = Specie("Cu")
        t = Specie("Ti")
        sites   = [[[c, c],
                    [t, t]],

                   [[c, t],
                    [t, c]]]
        sg = SitesGrid(sites)
        gen = CStru.gen_speckle(self.m, Specie("Cu"), (2,2,2), Specie("Ti"), 4)
        from collections import Iterator
        self.assertIsInstance(gen, Iterator)
        self.assertIn(CStru(self.m, sg), gen)
        self.assertEqual(next(gen).get_array().sum(), 204)
        self.assertEqual(next(gen).get_array().sum(), 204)

    def test_get_cell(self):
        c = Specie("Cu")
        t = Specie("Ti")
        m = [[-0.5, -0.5, -0.5],
             [-0.5,  0.5,  0.5],
             [ 0.5, -0.5,  0.5]]

        sites01 = [[[c]]]
        sites02 = [[[t, c, t],
                    [t, t, c]]]
        sg01 = SitesGrid(sites01)
        sg02 = SitesGrid(sites02)
        cstru01 = CStru(m, sg01)
        cstru02 = CStru(m, sg02)
        lat01, pos01, num01 = cstru01.get_cell()
        lat02, pos02, num02 = cstru02.get_cell()
        self.assertTrue(np.allclose(lat01, np.array([[-0.5, -0.5, -0.5],
                                                    [-0.5,  0.5,  0.5],
                                                    [ 0.5, -0.5,  0.5]])))
        self.assertTrue(np.allclose(pos01, np.array([[0, 0, 0]])))
        self.assertTrue(np.allclose(num01, np.array([29])))

        self.assertTrue(np.allclose(lat02, np.array([[-0.5, -0.5, -0.5],
                                          [-1,  1,  1],
                                         [ 1.5, -1.5,  1.5]])))
        self.assertTrue(np.allclose(pos02, np.array([[0, 0, 0],
                                         [0, 0, 1/3],
                                         [0, 0, 2/3],
                                         [0, 1/2, 0],
                                         [0, 1/2, 1/3],
                                         [0, 1/2, 2/3]])))
        self.assertTrue(np.allclose(num02, np.array([22, 29, 22, 22, 22, 29])))


class testGeneralCell(unittest.TestCase):

    def setUp(self):
        arr_lat = np.array([[3.0, 0, 0], [0, 2.0, 0.0], [0, 0, 1.0]])
        positions = [
                        [0.00000, 0.00000, 0.00000],
                        [0.00000, 0.50000, 0.00000],
                        [0.33333, 0.00000, 0.00000],
                        [0.33333, 0.50000, 0.00000],
                        [0.66666, 0.00000, 0.00000],
                        [0.66666, 0.50000, 0.00000],
                        [0.16666, 0.25000, 0.50000],
                        [0.16666, 0.75000, 0.50000],
                        [0.50000, 0.25000, 0.50000],
                        [0.50000, 0.75000, 0.50000],
                        [0.83333, 0.25000, 0.50000],
                        [0.83333, 0.75000, 0.50000]
                    ]
        arr_positions = np.array(positions)
        arr_numbers = np.array([6]*12)
        self.cell = GeneralCell(arr_lat, arr_positions, arr_numbers)

    # def test_eq(self):
    #     arr_lat = np.array([[3.0, 0, 0], [0, 2.0, 0.0], [0, 0, 1.0]])
    #     positions = [
    #                     [0.00000, 0.00000, 0.00000],
    #                     [0.00000, 0.50000, 0.00000],
    #                     [0.33333, 0.00000, 0.00000],
    #                     [0.33333, 0.50000, 0.00000],
    #                     [0.66666, 0.00000, 0.00000],
    #                     [0.66666, 0.50000, 0.00000],
    #                     [0.16666, 0.25000, 0.50000],
    #                     [0.16666, 0.75000, 0.50000],
    #                     [0.50000, 0.25000, 0.50000],
    #                     [0.50000, 0.75000, 0.50000],
    #                     [0.83333, 0.25000, 0.50000],
    #                     [0.83333, 0.75000, 0.50000]
    #                 ]
    #     arr_positions = np.array(positions)
    #     arr_numbers = np.array([6]*12)
    #     cell = GeneralCell(arr_lat, arr_positions, arr_numbers)
    #     self.assertEqual(cell, self.cell)
    #
    #     arr_lat = np.array([[3.0, 0, 0], [0, 2.0, 0.0], [0, 0, 1.0]])
    #     positions = [
    #                     [0.00000, 0.00000, 0.00000],
    #                     [0.00000, 0.50000, 0.00000],
    #                     [0.33333, 0.00000, 0.00000],
    #                     [0.66666, 0.50000, 0.00000],
    #                     [0.16666, 0.25000, 0.50000],
    #                     [0.16666, 0.75000, 0.50000],
    #                     [0.33333, 0.50000, 0.00000],
    #                     [0.66666, 0.00000, 0.00000],
    #                     [0.50000, 0.25000, 0.50000],
    #                     [0.50000, 0.75000, 0.50000],
    #                     [0.83333, 0.25000, 0.50000],
    #                     [0.83333, 0.75000, 0.50000]
    #                 ]
    #     arr_positions = np.array(positions)
    #     arr_numbers = np.array([6]*12)
    #     cell = GeneralCell(arr_lat, arr_positions, arr_numbers)
    #     self.assertEqual(cell, self.cell)
    #
    #     arr_lat = np.array([[3.0, 0, 0], [0, 2.0, 0.0], [0, 0, 1.0]])
    #     positions = [
    #                     [0.00000, 0.00000, 0.00000],
    #                     [0.00000, 0.50000, 0.00000],
    #                     [0.33333, 0.00000, 0.00000],
    #                     [0.33333, 0.50000, 0.00000],
    #                     [0.66666, 0.00000, 0.00000],
    #                     [0.66666, 0.50000, 0.00000],
    #                     [0.16666, 0.25000, 0.50000],
    #                     [0.16666, 0.75000, 0.50000],
    #                     [0.50000, 0.25000, 0.50000],
    #                     [0.50000, 0.75000, 0.50000],
    #                     [0.83333, 0.25000, 0.50000],
    #                     [0.83333, 0.75000, 0.50000]
    #                 ]
    #     arr_positions = np.array(positions)
    #     arr_numbers = np.array([4]*12)
    #     cell = GeneralCell(arr_lat, arr_positions, arr_numbers)
    #     self.assertNotEqual(cell, self.cell)
    #
    #     arr_lat = np.array([[2.0, 0, 0], [0, 2.0, 0.0], [0, 0, 1.0]])
    #     positions = [
    #                     [0.00000, 0.00000, 0.00000],
    #                     [0.00000, 0.50000, 0.00000],
    #                     [0.33333, 0.00000, 0.00000],
    #                     [0.33333, 0.50000, 0.00000],
    #                     [0.66666, 0.00000, 0.00000],
    #                     [0.66666, 0.50000, 0.00000],
    #                     [0.16666, 0.25000, 0.50000],
    #                     [0.16666, 0.75000, 0.50000],
    #                     [0.50000, 0.25000, 0.50000],
    #                     [0.50000, 0.75000, 0.50000],
    #                     [0.83333, 0.25000, 0.50000],
    #                     [0.83333, 0.75000, 0.50000]
    #                 ]
    #     arr_positions = np.array(positions)
    #     arr_numbers = np.array([6]*12)
    #     cell = GeneralCell(arr_lat, arr_positions, arr_numbers)
    #     self.assertNotEqual(cell, self.cell)

    def test_get_speckle_num(self):
        self.assertEqual(self.cell.get_speckle_num(Specie("B")), 0)
        self.assertEqual(self.cell.get_speckle_num(Specie("C")), 12)

    def test_property(self):
        np.testing.assert_equal(self.cell.numbers, np.array([6]*12))
        self.assertEqual(self.cell.comment, 'C12')


    def test_get_symmetry(self):
        sym = self.cell.get_spacegroup()
        self.assertEqual(sym, 'Im-3m (229)')

    def test_get_symmetry_permutation(self):
        sym_num = len(self.cell.get_symmetry_permutation())
        self.assertEqual(sym_num, 96)

    def test_is_primitive(self):
        self.assertFalse(self.cell.is_primitive())

    def test_degeneracy(self):
        sym_perm = self.cell.get_symmetry_permutation()
        self.assertEqual(self.cell.get_degeneracy(sym_perm), 1)

        ### more complex test
        arr_lat = np.array([[3.0, 0, 0], [0, 2.0, 0.0], [0, 0, 1.0]])
        positions = [
                        [0.00000, 0.00000, 0.00000],
                        [0.00000, 0.50000, 0.00000],
                        [0.33333, 0.00000, 0.00000],
                        [0.33333, 0.50000, 0.00000],
                        [0.66666, 0.00000, 0.00000],
                        [0.66666, 0.50000, 0.00000],
                        [0.16666, 0.25000, 0.50000],
                        [0.16666, 0.75000, 0.50000],
                        [0.50000, 0.25000, 0.50000],
                        [0.50000, 0.75000, 0.50000],
                        [0.83333, 0.25000, 0.50000],
                        [0.83333, 0.75000, 0.50000]
                    ]
        arr_positions = np.array(positions)
        arr_numbers = np.array([5,6,6,6,6,6,6,6,6,6,6,6])
        gcell = GeneralCell(arr_lat, arr_positions, arr_numbers)
        self.assertEqual(gcell.get_degeneracy(sym_perm), 12)

        arr_numbers = np.array([5,5,6,6,6,6,6,6,6,6,6,6])
        gcell = GeneralCell(arr_lat, arr_positions, arr_numbers)
        self.assertEqual(gcell.get_degeneracy(sym_perm), 6)

    def test_get_cartesian(self):
        arr_lat = np.array([[3.0, 0, 0], [0, 2.0, 0.0], [0, 0, 1.0]])
        positions = [
                        [0.00000, 0.00000, 0.00000],
                        [0.00000, 0.50000, 0.00000],
                        [0.33333, 0.00000, 0.00000],
                        [0.33333, 0.50000, 0.00000],
                        [0.66666, 0.00000, 0.00000],
                        [0.66666, 0.50000, 0.00000],
                        [0.16666, 0.25000, 0.50000],
                        [0.16666, 0.75000, 0.50000],
                        [0.50000, 0.25000, 0.50000],
                        [0.50000, 0.75000, 0.50000],
                        [0.83333, 0.25000, 0.50000],
                        [0.83333, 0.75000, 0.50000]
                    ]
        arr_positions = np.array(positions)
        arr_numbers = np.array([6]*12)
        cell = GeneralCell(arr_lat, arr_positions, arr_numbers)

        cart_coor = cell.get_cartesian()
        cart_ans = np.array([[ 0.     ,  0.     ,  0.     ],
                             [ 0.     ,  1.     ,  0.     ],
                             [ 0.49998,  0.5    ,  0.5    ],
                             [ 0.49998,  1.5    ,  0.5    ],
                             [ 0.99999,  0.     ,  0.     ],
                             [ 0.99999,  1.     ,  0.     ],
                             [ 1.5    ,  0.5    ,  0.5    ],
                             [ 1.5    ,  1.5    ,  0.5    ],
                             [ 1.99998,  0.     ,  0.     ],
                             [ 1.99998,  1.     ,  0.     ],
                             [ 2.49999,  0.5    ,  0.5    ],
                             [ 2.49999,  1.5    ,  0.5    ]])
        self.assertTrue(np.allclose(cart_coor, cart_ans))

        arr_numbers = np.array([6, 6, 6, 5, 6, 5, 6, 6, 6, 6, 6, 6])
        cell = GeneralCell(arr_lat, arr_positions, arr_numbers)

        cart_coor = cell.get_cartesian(Specie('B'))
        cart_ans = np.array([[ 0.49998,  1.5    ,  0.5    ],
                             [ 0.99999,  1.     ,  0.     ]])



    def test_supercell(self):
        scale_mat = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
        scell = self.cell.supercell(scale_mat)
        arr_lat = np.array([[3.0, 0, 0], [0, 2.0, 0.0], [0, 0, 2.0]])
        positions =  [[ 0.     ,  0.     ,  0.     ],
                      [ 0.     ,  0.     ,  0.5    ],
                      [ 0.     ,  0.5    ,  0.     ],
                      [ 0.     ,  0.5    ,  0.5    ],
                      [ 0.16666,  0.25   ,  0.25   ],
                      [ 0.16666,  0.25   ,  0.75   ],
                      [ 0.16666,  0.75   ,  0.25   ],
                      [ 0.16666,  0.75   ,  0.75   ],
                      [ 0.33333,  0.     ,  0.     ],
                      [ 0.33333,  0.     ,  0.5    ],
                      [ 0.33333,  0.5    ,  0.     ],
                      [ 0.33333,  0.5    ,  0.5    ],
                      [ 0.5    ,  0.25   ,  0.25   ],
                      [ 0.5    ,  0.25   ,  0.75   ],
                      [ 0.5    ,  0.75   ,  0.25   ],
                      [ 0.5    ,  0.75   ,  0.75   ],
                      [ 0.66666,  0.     ,  0.     ],
                      [ 0.66666,  0.     ,  0.5    ],
                      [ 0.66666,  0.5    ,  0.     ],
                      [ 0.66666,  0.5    ,  0.5    ],
                      [ 0.83333,  0.25   ,  0.25   ],
                      [ 0.83333,  0.25   ,  0.75   ],
                      [ 0.83333,  0.75   ,  0.25   ],
                      [ 0.83333,  0.75   ,  0.75   ]]
        arr_positions = np.array(positions)
        arr_numbers = np.array([6]*24)
        ans_scell = GeneralCell(arr_lat, arr_positions, arr_numbers)
        self.assertTrue(np.allclose(scell.positions, arr_positions))


class testModifiedCell(unittest.TestCase):


    def setUp(self):
        self.latt = np.array([[12.,2.,14.],[16.,0.,16.],[2.,2.,0.]])
        self.pos = np.array([[0.0, 0.0, 0.0],
                             [0.0, 0.125, 0.0],
                             [0.0, 0.25, 0.0],
                             [0.0, 0.375, 0.0],
                             [0.0, 0.5, 0.0],
                             [0.0, 0.625, 0.0],
                             [0.0, 0.75, 0.0],
                             [0.0, 0.875, 0.0],
                             [0.25, 0.09375, 0.25],
                             [0.25, 0.21875, 0.25],
                             [0.25, 0.34375, 0.25],
                             [0.25, 0.46875, 0.25],
                             [0.25, 0.59375, 0.25],
                             [0.25, 0.71875, 0.25],
                             [0.25, -0.15625, 0.25],
                             [0.25, -0.03125, 0.25]]
                             )
        self.numbers = np.array([30, 30, 30, 30, 30, 30, 30, 30, 16, 16, 16, 16, 16, 16, 16, 16])
        self.initstru = ModifiedCell(self.latt)
        self.full_initzb = ModifiedCell(self.latt, self.pos, self.numbers)

    def test_from_to_gcell(self):
        gcell = GeneralCell(self.latt, self.pos, self.numbers)
        modcell = ModifiedCell.from_gcell(gcell)
        new_gcell = modcell.to_gcell()
        self.assertTrue(np.allclose(new_gcell.positions, gcell.positions))

    def test_get_cartesian_from_frac(self):
        expect_coor = np.array([[  0.,   0.,   0.],
                                 [  2.,   0.,   2.],
                                 [  4.,   0.,   4.],
                                 [  6.,   0.,   6.],
                                 [  8.,   0.,   8.],
                                 [ 10.,   0.,  10.],
                                 [ 12.,   0.,  12.],
                                 [ 14.,   0.,  14.],
                                 [  5.,   1.,   5.],
                                 [  7.,   1.,   7.],
                                 [  9.,   1.,   9.],
                                 [ 11.,   1.,  11.],
                                 [ 13.,   1.,  13.],
                                 [ 15.,   1.,  15.],
                                 [  1.,   1.,   1.],
                                 [  3.,   1.,   3.]])

        cart_coor = self.full_initzb.get_cartesian_from_frac(self.pos)
        self.assertTrue(np.allclose(cart_coor, expect_coor))

    def test_get_frac_from_cart(self):
        """
        change a to b and then b to a
        make sure that a equals a
        """
        expect_coor = self.pos
        cart_coor = self.full_initzb.get_cartesian_from_frac(self.pos)
        frac_coor = self.full_initzb.get_frac_from_cart(cart_coor)
        self.assertTrue(np.allclose(frac_coor, expect_coor))

    def test_translate_sites(self):
        def get_rand_vec():
            # deals with zero vectors.
            vector = np.random.randn(3)
            vnorm = np.linalg.norm(vector)
            return vector / vnorm * 0.2 if vnorm != 0 else get_rand_vec()

        stru = self.full_initzb.copy()
        stru.translate_sites(0, get_rand_vec())
        gcell = stru.to_gcell()
        # This is a runabal passive test TODO!!!!
        self.assertFalse(np.array_equal(gcell.positions, self.pos))

    def test_perturb(self):
        stru = self.full_initzb.copy()
        stru.perturb(0.1)
        gcell = stru.to_gcell()
        # This is a runabal passive test TODO!!!!
        self.assertFalse(np.array_equal(gcell.positions, self.pos))

        stru = self.full_initzb.copy()
        stru.perturb(0.0)
        gcell = stru.to_gcell()
        # This is a runabal passive test TODO!!!!
        self.assertTrue(np.array_equal(gcell.positions, self.pos))

    def test_swap_axis(self):
        lattice = np.array([[3.464, 3.0, 0.0],
                           [2.0, 5.99, 0.0],
                           [0.0, 0.0, 15]])
        positions = np.array([[0.014542,  0.690903,   0.508516],
                              [0.019469,  0.360318,   0.506407],
                              [0.004450,  0.027941,   0.507198],
                              [0.464436,  0.525638,   0.502271],
                              [0.523634,  0.136503,   0.501554],
                              [0.540583,  0.822100,   0.508352]])
        numbers = np.array([4,4,6,6,6,6])
        modcell = ModifiedCell(lattice, positions, numbers)
        modcell.swap_axis((0,2,1))
        expect_lattice = np.array([[3.464, 0.0, 3.0],
                                   [0.0, 15, 0.0],
                                   [2.0, 0.0, 5.99]])
        expect_positions = np.array([[0.014542, 0.508516, 0.690903],
                                     [0.019469, 0.506407, 0.360318],
                                     [0.004450, 0.507198, 0.027941],
                                     [0.464436, 0.502271, 0.525638],
                                     [0.523634, 0.501554, 0.136503],
                                     [0.540583, 0.508352, 0.822100]])
        expect_numbers = np.array([4,4,6,6,6,6])
        self.assertTrue(np.array_equal(modcell.lattice, expect_lattice))
        self.assertTrue(np.array_equal(modcell.positions, expect_positions))
        self.assertTrue(np.array_equal(modcell.numbers, expect_numbers))

        # lattice = np.array([[3.464, 3.0, 0.0],
        #                    [2.0, 5.99, 0.0],
        #                    [0.0, 0.0, 15]])
        # positions = np.array([[0.014542,  0.690903,   0.508516],
        #                       [0.019469,  0.360318,   0.506407],
        #                       [0.004450,  0.027941,   0.507198],
        #                       [0.464436,  0.525638,   0.502271],
        #                       [0.523634,  0.136503,   0.501554],
        #                       [0.540583,  0.822100,   0.508352]])
        # numbers = np.array([4,4,6,6,6,6])
        # modcell = ModifiedCell(lattice, positions, numbers)
        # modcell.swap_axis((2,0,1))
        # expect_lattice = np.array([[0.0, 0.0, 15],
        #                            [3.464, 3.0, 0.0],
        #                            [2.0, 5.99, 0.0]])
        # expect_positions = np.array([[0.508516, 0.014542, 0.690903],
        #                              [0.506407, 0.019469, 0.360318],
        #                              [0.507198, 0.004450, 0.027941],
        #                              [0.502271, 0.464436, 0.525638],
        #                              [0.501554, 0.523634, 0.136503],
        #                              [0.508352, 0.540583, 0.822100]])
        # expect_numbers = np.array([4,4,6,6,6,6])
        # self.assertTrue(np.array_equal(modcell.lattice, expect_lattice))
        # self.assertTrue(np.array_equal(modcell.positions, expect_positions))
        # self.assertTrue(np.array_equal(modcell.numbers, expect_numbers))

    def test_d2_at_Z(self):
        lattice = np.array([[5.29, 0.0, 0.0],
                           [0.0, 15.0, 0.0],
                           [-0.756071, 0.0, 5.237142]])
        positions = np.array([[0.286607, 0.490220, 0.435572],
                              [0.414891, 0.496366, 0.780603],
                              [0.546188, 0.504194, 0.217210],
                              [0.681347, 0.501510, 0.598678],
                              [0.090962, 0.507718, 0.704071],
                              [0.182893, 0.512417, 0.074836],
                              [0.823989, 0.511788, 0.925624],
                              [0.899757, 0.501391, 0.302911]])
        numbers = np.array([4,4,4,4,6,6,6,6])
        modcell = ModifiedCell(lattice, positions, numbers)
        modcell.d2_at_Z(z=15.0)
        expect_lattice = np.array([[5.29, 0.0, 0.0],
                                   [-0.756071, 5.237142, 0.0],
                                   [0.0, 0.0, 15.0]])
        expect_positions = np.array([[0.286607, 0.435572, 0.490220],
                                     [0.414891, 0.780603, 0.496366],
                                     [0.546188, 0.217210, 0.504194],
                                     [0.681347, 0.598678, 0.501510],
                                     [0.090962, 0.704071, 0.507718],
                                     [0.182893, 0.074836, 0.512417],
                                     [0.823989, 0.925624, 0.511788],
                                     [0.899757, 0.302911, 0.501391]])
        expect_numbers = np.array([4,4,4,4,6,6,6,6])
        self.assertTrue(np.array_equal(modcell.lattice, expect_lattice))
        self.assertTrue(np.array_equal(modcell.positions, expect_positions))
        self.assertTrue(np.array_equal(modcell.numbers, expect_numbers))

    def test_get_points_in_sphere(self):
        latt = np.array([[4.898979, 0.000000, 0.000000],
                         [2.449490,  4.242641, 0.000000],
                         [1.632993, -0.000000, 4.618802]])
        pos = np.array([[0.208333, 0.333333, 0.375000],
                        [0.375000, 0.000000, 0.875000],
                        [0.541667, 0.666667, 0.375000],
                        [0.708333, 0.333333, 0.875000],
                        [0.875000, 0.000000, 0.375000],
                        [0.000000, 0.000000, 0.000000],
                        [0.166667, 0.666667, 0.500000],
                        [0.333333, 0.333333, 0.000000],
                        [0.500000, 0.000000, 0.500000],
                        [0.666667, 0.666667, 0.000000],
                        [0.833333, 0.333333, 0.500000],
                        [0.041667, 0.666667, 0.875000]])
        numbers = np.array([16,16,16,16,16,30,30,30,30,30,30,55])
        modcell = ModifiedCell(latt, pos, numbers)
        dict_sites = modcell.get_points_incell_insphere(np.array([0.041667, 0.666667, 0.875000]), 2)

        sites = [Site(pos[5], 'Zn'),
                 Site(pos[9], 'Zn'),
                 Site(pos[7], 'Zn'),
                 Site(pos[6], 'Zn'),
                 Site(pos[11], 'Cs')]
        self.assertEqual(len(dict_sites), 5)
        for s in dict_sites.values():
            self.assertTrue(s in sites)

        #  Test for find giving element
        modcell = ModifiedCell(latt, pos, numbers)
        dict_sites = modcell.get_points_incell_insphere(np.array([0.041667, 0.666667, 0.875000]), 2, Specie('Zn'))

        sites = [Site(pos[5], 'Zn'),
                 Site(pos[9], 'Zn'),
                 Site(pos[7], 'Zn'),
                 Site(pos[6], 'Zn')]
        self.assertEqual(len(dict_sites), 4)
        for s in dict_sites.values():
            self.assertTrue(s in sites)

    def test_append(self):
        site_s = Site((0,0,0.222), 'S')
        zb = self.full_initzb.copy()
        zb.append(site_s)
        gcell = zb.to_gcell()
        new_pos = np.insert(self.pos, 1, [[0,0,0.222]], axis=0)
        new_numbers = np.insert(self.numbers, 1, 16)
        self.assertTrue(np.array_equal(gcell.positions, new_pos))
        self.assertTrue(np.array_equal(gcell.numbers, new_numbers))

    def test_remove_site(self):
        stru = self.full_initzb.copy()
        new_stru = stru.remove_site()
        gcell = new_stru.to_gcell()
        new_pos = np.delete(self.pos, 15, 0)
        self.assertTrue(np.array_equal(gcell.positions, new_pos))

        stru = self.full_initzb.copy()
        new_stru = stru.remove_site(2)
        gcell = new_stru.to_gcell()
        new_pos = np.delete(self.pos, 2, 0)
        self.assertTrue(np.array_equal(gcell.positions, new_pos))

    def test_append_site(self):
        site_s = Site((0,0,0.222), 'S')
        initstru = self.initstru.copy()
        new_stru = initstru.append_site(site_s)
        gcell = new_stru.to_gcell()
        new_pos = np.append([[0,0,0]], [[0,0,0.222]], axis=0)
        self.assertTrue(new_stru is initstru)
        self.assertTrue(np.array_equal(gcell.positions, new_pos))

    def test_remove_sites(self):
        stru = self.full_initzb.copy()
        new_stru = stru.remove_sites([2, 14, 3, 5])
        gcell = new_stru.to_gcell()
        new_pos = np.delete(self.pos, [2, 14, 3, 5], 0)
        self.assertTrue(np.array_equal(gcell.positions, new_pos))

    def test_append_sites(self):
        l_sites = [Site((0,0,k), 'S') for k in [0.1, 0.2, 0.3]]
        initstru = self.initstru.copy()
        new_stru = initstru.append_sites(l_sites)
        gcell = new_stru.to_gcell()
        new_pos = np.array([[0,0,0],[0,0,0.1],[0,0,0.2],[0,0,0.3]])
        self.assertTrue(np.array_equal(gcell.positions, new_pos))

    def test_copy(self):
        new_zb = self.full_initzb.copy()
        self.assertEqual(new_zb, self.full_initzb)
        self.assertFalse(new_zb is self.full_initzb)


if __name__ == "__main__":
    import nose2
    nose2.main()
