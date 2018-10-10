# -*- coding: utf-8 -*-
import unittest
import numpy

from sagar.crystal.structure import Cell
from sagar.crystal.derive import ConfigurationGenerator as CG


class TestDerive(unittest.TestCase):

    def setUp(self):
        fcc_latt = [0, 5, 5,
                    5, 0, 5,
                    5, 5, 0]
        fcc_pos = [(0, 0, 0), (0.5, 0.5, 0.5)]
        fcc_atoms = [1, 2]
        self.fcc_pcell = Cell(fcc_latt, fcc_pos, fcc_atoms)

    def test_cons_max_volume(self):
        # 最大体积下可能产生的所有结构的总和，不包含超胞。
        wanted = [2, 4, 10, 29]
        got = []
        cg = CG(self.fcc_pcell)
        for v in [1, 2, 3, 4]:
            con = cg.cons_max_volume([(1, 5), (2,)], v)
            got.append(len([i for i in con]))

        self.assertEqual(got, wanted)

    def test_cons_specific_volume(self):
        wanted = [2, 6, 12, 41]
        got = []
        cg = CG(self.fcc_pcell)
        for v in [1, 2, 3, 4]:
            con = cg.cons_specific_volume([(1, 5), (2,)], v)
            got.append(len([i for i in con]))

        self.assertEqual(got, wanted)

    def test_cons_specific_volume_c(self):
        wanted = [7]
        got = []
        cg = CG(self.fcc_pcell)
        for v in [4]:
            con = cg.cons_specific_volume([(1, 5), (2,)], v, e_num=(3, 1))
            got.append(len([i for i in con]))

        self.assertEqual(got, wanted)

    # def test_degeneracy_of_confs_nondup_specific_volume(self):
    #     wanted = [1, 2, 4, 16, 32]
    #     got = []
    #     for v in [1, 2, 3, 4, 5]:
    #         deg_all = 0
    #         cons = confs_nondup_specific_volume(self.fcc_pcell, [[1, 5], [2]], v)
    #         for c in cons:
    #             deg_all += c[1]
    #         print(deg_all)

    def test_cons_specific_cell(self):
        fcc_latt = [5, 0, 0,
                    0, 5, 0,
                    0, 0, 5]
        fcc_pos = [(0, 0, 0),
                   (0, 0.5, 0.5),
                   (0.5, 0, 0.5),
                   (0.5, 0.5, 0)]
        fcc_atoms = [0, 0, 0, 0]
        con_cell = Cell(fcc_latt, fcc_pos, fcc_atoms)
        cg = CG(con_cell)
        con = cg.cons_specific_cell([(2, 3), (2, 3), (2, 3), (2, 3)])

        # number of all configurations
        wanted = 5
        got = len([i for i in con])

        self.assertEqual(got, wanted)

    def test_cons_specific_cell_and_c(self):
        # FCC conventional cell
        fcc_latt = [5, 0, 0,
                    0, 5, 0,
                    0, 0, 5]
        fcc_pos = [(0, 0, 0),
                   (0, 0.5, 0.5),
                   (0.5, 0, 0.5),
                   (0.5, 0.5, 0)]
        fcc_atoms = [0, 0, 0, 0]
        con_cell = Cell(fcc_latt, fcc_pos, fcc_atoms)
        cg = CG(con_cell)
        con = cg.cons_specific_cell(
            [(2, 3, 4), (2, 3, 4), (2, 3, 4), (2, 3, 4)], e_num=(2, 1, 1))

        # number of all configurations
        wanted = 1
        got = len([i for i in con])

        self.assertEqual(got, wanted)

        # Zinc-blende conventional cell
        fcc_latt = [5, 0, 0,
                    0, 5, 0,
                    0, 0, 5]
        fcc_pos = [(0, 0, 0),
                   (0, 0.5, 0.5),
                   (0.5, 0, 0.5),
                   (0.5, 0.5, 0),
                   (0.25, 0.25, 0.25),
                   (0.75, 0.75, 0.25),
                   (0.75, 0.25, 0.75),
                   (0.25, 0.75, 0.75)]
        fcc_atoms = [0, 0, 0, 0, 1, 1, 1, 1]
        con_cell = Cell(fcc_latt, fcc_pos, fcc_atoms)
        cg = CG(con_cell)
        con = cg.cons_specific_cell(
            [(2, 3, 4), (2, 3, 4), (2, 3, 4), (2, 3, 4), (1, ), (1, ), (1, ), (1,)], e_num=(2, 1, 1))

        # number of all configurations
        wanted = 1
        got = len([i for i in con])

        self.assertEqual(got, wanted)


if __name__ == "__main__":
    import nose2
    nose2.main()
