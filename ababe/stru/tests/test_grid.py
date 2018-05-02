# coding: utf-8
# Distributed under the terms of the MIT License.
#
import unittest
import numpy as np

from ababe.stru.grid import SuperLatticeGenerator, SuperLatticeCell
from ababe.stru.grid import SuperLatticeGenerator2D
from ababe.stru.scaffold import GeneralCell


class testSuperLatticeGenerator(unittest.TestCase):

    def setUp(self):
        self.bcc_base = np.array([[0.5, 0.5, -0.5],
                                  [-0.5, 0.5, 0.5],
                                  [0.5, -0.5, 0.5]])
        self.bcc_u_position = np.array([[0, 0, 0]])
        self.bcc_u_n = np.array([0])
        self.bcc_uc = GeneralCell(self.bcc_base,
                                  self.bcc_u_position,
                                  self.bcc_u_n)

    def test_HNFs_from_n_dups(self):
        hnfs = SuperLatticeGenerator.hnfs_from_n_dups(self.bcc_uc, 4)
        self.assertEqual(len(hnfs), 35)

        hnfs = SuperLatticeGenerator.hnfs_from_n_dups(self.bcc_uc, 5)
        self.assertEqual(len(hnfs), 31)

        hnfs = SuperLatticeGenerator.hnfs_from_n_dups(self.bcc_uc, 6)
        self.assertEqual(len(hnfs), 91)

    def test_eq(self):
        hnfs = SuperLatticeGenerator.hnfs_from_n_dups(self.bcc_uc, 2)
        hnf00 = hnfs[0]
        hnf01 = hnfs[1]
        self.assertEqual(hnf00, hnf01)

    def test_HNFs_from_n(self):
        # use bcc data from <PHYSICAL REVIEW B 80, 014120 (2009)>
        # to test...
        results = [2, 3, 7, 5, 10, 7]
        for i, result in zip(range(2, 8), results):
            nodup_hnfs = SuperLatticeGenerator.hnfs_from_n(self.bcc_uc, i)
            self.assertEqual(len(nodup_hnfs), result)

        # hcp test
        hcp_b = np.array([[2.51900005,  0.,  0.],
                          [-1.25950003,  2.18151804, 0.],
                          [0., 0.,  4.09100008]])
        hcp_positions = np.array([[0.33333334,  0.66666669,  0.25],
                                  [0.66666663,  0.33333331,  0.75]])
        hcp_numbers = np.array([0, 0])
        hcp_uc = GeneralCell(hcp_b, hcp_positions, hcp_numbers)
        # import pdb; pdb.set_trace()
        result = [len(SuperLatticeGenerator.hnfs_from_n(hcp_uc, i)) for i in range(1,10)]
        self.assertEqual(result, [1, 3, 5, 11, 7, 19, 11, 34, 23])

    def test_to_general_cell(self):
        hnfs = SuperLatticeGenerator.hnfs_from_n_dups(self.bcc_uc, 4)
        hnf = hnfs[2]
        cell = hnf.to_general_cell()
        # print(cell.spg_cell)
        self.assertEqual(len(cell.spg_cell[2]), 4)

        zb_b = np.array([[3.82863, 0., 0.],
                         [1.91431, 3.31569, 0.],
                         [1.91431, 1.10523, 3.12606]])
        zb_pos = np.array([[0., 0., 0.],
                           [0.25, 0.25, 0.25]])
        zb_num = np.array([30, 16])
        zb_uc = GeneralCell(zb_b, zb_pos, zb_num)
        zbs = SuperLatticeGenerator.hnfs_from_n(zb_uc, 3)
        one_zb = zbs[1]
        cell = one_zb.to_general_cell()
        # print(cell.spg_cell)
        self.assertEqual(len(cell.spg_cell[2]), 6)

        h = np.array([[1, 0, 3],
                      [0, 2, 2],
                      [0, 0, 4]])
        strange_hnf = SuperLatticeCell(self.bcc_uc, h)
        shnf = strange_hnf.to_general_cell()
        self.assertEqual(len(shnf.spg_cell[2]), 8)


class testSuperLatticeGenerator2D(unittest.TestCase):

    def test_test(self):
        zb_b = np.array([[3.82863, 0., 0.],
                         [1.91431, 3.31569, 0.],
                         [1.91431, 1.10523, 3.12606]])
        zb_pos = np.array([[0., 0., 0.],
                           [0.25, 0.25, 0.25]])
        zb_num = np.array([30, 16])
        zb_uc = GeneralCell(zb_b, zb_pos, zb_num)
        # s2d = SuperLatticeGenerator2D(zb_uc, 1)

    def test_fatal(self):
        """
        make sure the number of output is right
        """
        pass



if __name__ == "__main__":
    import nose2
    nose2.main()
