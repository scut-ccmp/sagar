import unittest
import numpy

from pyabc.crystal.structure import Cell, MutableCell, frac_to_car
from pyabc.crystal.mutate import *
from pyabc.utils.math import distance


class TestMutableMethods(unittest.TestCase):

    def setUp(self):
        # 使用金刚石Si作为测试模型
        self.lattice = numpy.array([0.0, 0.5, 0.5,
                                    0.5, 0.0, 0.5,
                                    0.5, 0.5, 0.0]).reshape((3, 3))

    def test_cell_to_mcell(self):
        si_pos = [(-0.125, -0.125, -0.125),
                   (0.125, 0.125, 0.125)]
        si_atoms = [14, 14]
        c = Cell(self.lattice, si_pos, si_atoms)
        mc = cell_to_mcell(c)
        self.assertTrue(isinstance(mc, MutableCell))

        # 转换回去，并测试返回后的数据结构
        c = mc.to_cell()
        numpy.testing.assert_almost_equal(c.positions, numpy.array([0.875, 0.875, 0.875,
                                                                    0.125, 0.125, 0.125]).reshape((2, 3)))
        numpy.testing.assert_almost_equal(c.atoms, numpy.array([14, 14]))

    def test_perturb(self):
        lattice = numpy.copy(self.lattice)
        si_sites = [[(-0.125, -0.125, -0.125), "Si"],
                    [(0.125, 0.125, 0.125), "Si"]]
        import copy
        og_si_sites = copy.deepcopy(si_sites)
        mcell = MutableCell(lattice, sites=si_sites)
        perturb(mcell, distance=0.01)
        self.assertEqual(len(mcell._sites), 2)

        for i in range(2):
            og_frac_pos = numpy.array(og_si_sites[i][0])
            og_car_pos = frac_to_car(lattice, og_frac_pos)

            new_frac_pos = numpy.array(mcell._sites[i][0])
            new_car_pos = frac_to_car(lattice, new_frac_pos)

            numpy.testing.assert_almost_equal(distance(og_car_pos, new_car_pos), 0.01)

    def test_remove_sites_in_a_circle(self):
        # 需要测试两类情况，范围内只包含胞内原子的，和范围超出一个胞。
        # 只包含胞内原子

        # 范围超出一个胞，圆心在晶胞六面体顶点的特殊情况
        pass
