# -*- coding: utf-8 -*-
import unittest
import numpy
import copy

from sagar.crystal.structure import Cell, MutableCell, frac_to_car
from sagar.crystal.mutate import *
from sagar.toolkit.mathtool import distance


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

        # TODO: 需要确认产生的mcell的改变是不会影响原有的cell的

    def test_perturb(self):
        lattice = numpy.copy(self.lattice)
        si_sites = [[(-0.125, -0.125, -0.125), "Si"],
                    [(0.125, 0.125, 0.125), "Si"]]

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

    def test_rotate_site_by_z(self):
        lattice = numpy.array([4.912, 0.000, 0.000,
                               -2.456, 4.254, 0.000,
                               0.000, 0.000, 16.000])
        si_sites = [[(0.0, 0.0, 0.5), "C"],
                    [(0.0, 0.5, 0.5), "C"],
                    [(0.5, 0.0, 0.5), "C"],
                    [(0.5, 0.5, 0.5), "C"],
                    [(0.33333, 0.166666, 0.5), "C"],
                    [(0.33333, 0.666666, 0.5), "C"],
                    [(0.83333, 0.166666, 0.5), "C"],
                    [(0.83333, 0.666666, 0.5), "C"]]

        mcell = MutableCell(lattice, sites=si_sites)
        rotate_sites_in_a_circle_by_z(mcell, (0.416666, 0.5833333, 0.5), radius=1.0, degrees=90)

        c = mcell.to_cell()
        expected_positions = numpy.array([[0.       , 0.       , 0.5      ],
                                          [0.       , 0.5      , 0.5      ],
                                          [0.5      , 0.       , 0.5      ],
                                          [0.3333327, 0.4999993, 0.5      ],
                                          [0.33333  , 0.166666 , 0.5      ],
                                          [0.4999987, 0.6666693, 0.5      ],
                                          [0.83333  , 0.166666 , 0.5      ],
                                          [0.83333  , 0.666666 , 0.5      ]]
)
        numpy.testing.assert_almost_equal(c.positions, expected_positions)



    def test_remove_sites_in_a_circle(self):
        # 需要测试两类情况，范围内只包含胞内原子的，和范围超出一个胞。

        # 原始的晶体
        lattice = numpy.array([5.5, 0, 0,
                               0, 5.5, 0,
                               0, 0, 11.0]).reshape((3, 3))
        sites = [[(0.100, 0.100, 0.050), "Zn"],
                 [(0.100, 0.100, 0.550), "Zn"],
                 [(0.100, 0.600, 0.300), "Zn"],
                 [(0.100, 0.600, 0.800), "Zn"],
                 [(0.600, 0.100, 0.300), "Zn"],
                 [(0.600, 0.100, 0.800), "Zn"],
                 [(0.600, 0.600, 0.050), "Zn"],
                 [(0.600, 0.600, 0.550), "Zn"],
                 [(0.350, 0.350, 0.175), "S"],
                 [(0.350, 0.350, 0.675), "S"],
                 [(0.850, 0.850, 0.175), "S"],
                 [(0.850, 0.850, 0.675), "S"],
                 [(0.850, 0.350, 0.425), "S"],
                 [(0.850, 0.350, 0.925), "S"],
                 [(0.350, 0.850, 0.425), "S"],
                 [(0.350, 0.850, 0.925), "S"]]

        # Condition I
        mcell = MutableCell(numpy.copy(lattice), sites=copy.deepcopy(sites))

        # 只去除包含胞内原子
        # [(0.600, 0.600, 0.550), "Zn"] 该位点周围四个S原子都在胞内
        # Zn-S 键长2.38157A, 所以选择半径2.5A内作为截断半径
        remove_sites_in_a_circle(mcell, (0.6, 0.6, 0.55), 2.5, list_ele=['S'])

        # 预期去除以下S原子
        # [(0.350, 0.350, 0.675), "S"],
        # [(0.850, 0.850, 0.675), "S"],
        # [(0.850, 0.350, 0.425), "S"],
        # [(0.350, 0.850, 0.425), "S"],
        expected_sites = [[(0.100, 0.100, 0.050), "Zn"],
                          [(0.100, 0.100, 0.550), "Zn"],
                          [(0.100, 0.600, 0.300), "Zn"],
                          [(0.100, 0.600, 0.800), "Zn"],
                          [(0.600, 0.100, 0.300), "Zn"],
                          [(0.600, 0.100, 0.800), "Zn"],
                          [(0.600, 0.600, 0.050), "Zn"],
                          [(0.600, 0.600, 0.550), "Zn"],
                          [(0.350, 0.350, 0.175), "S"],
                          [(0.350, 0.350, 0.675), "Vacc"],
                          [(0.850, 0.850, 0.175), "S"],
                          [(0.850, 0.850, 0.675), "Vacc"],
                          [(0.850, 0.350, 0.425), "Vacc"],
                          [(0.850, 0.350, 0.925), "S"],
                          [(0.350, 0.850, 0.425), "Vacc"],
                          [(0.350, 0.850, 0.925), "S"]]
        self.assertEqual(mcell._sites, expected_sites)

        # Condition II
        mcell = MutableCell(numpy.copy(lattice), sites=copy.deepcopy(sites))
        # 范围超出一个胞，圆心在晶胞六面体顶点的特殊情况
        # 去除最靠近原点的Zn周围的S原子
        # [[(0.100, 0.100, 0.050), "Zn"],
        remove_sites_in_a_circle(mcell, (0.1, 0.1, 0.05), 2.5, list_ele=['S'])

        # 预期去除以下S原子
        # [(0.350, 0.350, 0.175), "S"],
        # [(0.850, 0.350, 0.925), "S"],
        # [(0.350, 0.850, 0.925), "S"],
        # [(0.850, 0.850, 0.175), "S"],
        expected_sites = [[(0.100, 0.100, 0.050), "Zn"],
                          [(0.100, 0.100, 0.550), "Zn"],
                          [(0.100, 0.600, 0.300), "Zn"],
                          [(0.100, 0.600, 0.800), "Zn"],
                          [(0.600, 0.100, 0.300), "Zn"],
                          [(0.600, 0.100, 0.800), "Zn"],
                          [(0.600, 0.600, 0.050), "Zn"],
                          [(0.600, 0.600, 0.550), "Zn"],
                          [(0.350, 0.350, 0.175), "Vacc"],
                          [(0.350, 0.350, 0.675), "S"],
                          [(0.850, 0.850, 0.175), "Vacc"],
                          [(0.850, 0.850, 0.675), "S"],
                          [(0.850, 0.350, 0.425), "S"],
                          [(0.850, 0.350, 0.925), "Vacc"],
                          [(0.350, 0.850, 0.425), "S"],
                          [(0.350, 0.850, 0.925), "Vacc"]]
        self.assertEqual(mcell._sites, expected_sites)
