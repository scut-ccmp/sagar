# -*- coding: utf-8 -*-

# 该模块包含针对`sagar.crystal.structure`中
# 的MutableCell产生的对象的操作。
import numpy
from itertools import product

from sagar.crystal.structure import MutableCell
from sagar.crystal.structure import get_symbol, car_to_frac, frac_to_car
from sagar.toolkit.mathtool import distance


def cell_to_mcell(cell):
    lattice = numpy.copy(cell.lattice)
    sites = []
    for pos, ele in zip(cell.positions.tolist(), cell.atoms.tolist()):
        symbol = get_symbol(ele)
        sites.append([pos] + [symbol])
    return MutableCell(cell.lattice, sites)


def perturb(mcell, distance=0.02):
    """
    对一个结构的所有原子进行扰动，扰动的距离是distance，单位是A
    distance默认为0.02A
    perturb sites one by one
    """
    for site in mcell._sites:
        vec = _get_rand_vec(distance)
        vec = car_to_frac(mcell._lattice, vec)
        site[0] = numpy.array(site[0]) + vec


def _get_rand_vec(distance):
    # deals with zero vectors.
    vector = numpy.random.randn(3)
    vnorm = numpy.linalg.norm(vector)
    return vector / vnorm * distance if vnorm != 0 else _get_rand_vec()


def remove_sites_in_a_circle(mcell, cc, radius, list_ele=None):
    """
    移除以cc为圆心，radius为半径的圆内，所有ele的原子。
    cc: 圆心 (x, y, z) 分数坐标
    radius: 半径 单位为A
    ele: 要删除的元素的list，若为None，则删除所有元素找到的位点
        这里使用元素符号的列表
    """
    for idx, s in enumerate(mcell._sites):
        if _is_close_in_radius(mcell._lattice, cc, s[0], radius) and _is_in_ele(s[1], list_ele):
            mcell.remove_site(idx)

def rotate_sites_in_a_circle_by_z(mcell, cc, radius, degrees, radians=None, list_ele=None):
    """
    cc: 圆心 (x, y, z) 分数坐标, 旋转时z不动
    radius: 半径范围内的原子 单位为A
    ele: 要旋转的元素的list，若为None，则旋转所有元素找到的位点
        这里使用元素符号的列表
    """
    for idx, s in enumerate(mcell._sites):
        if _is_close_in_radius(mcell._lattice, cc, s[0], radius) and _is_in_ele(s[1], list_ele):
            mcell.rotate_site_by_z(idx, cc, degrees)

def _is_close_in_radius(lattice, p1, p2, radius):
    # 截断圆对应的半径大于距离，说明两点在圆内，需要删除：返回True
    # 因为周期性边界条件，首先将原始胞沿着x,y,z三个方向得到27个位置坐标

    # TODO: 扩胞为3x3x3是不必要的，可以根据radius和lattice选择性扩胞，
    # 若lattice的边或顶点在半径以内才需要选择向该方向扩胞
    symprec = 1e-5

    trans = numpy.array([i for i in product([-1, 0, 1], repeat=3)])
    all_p1 = numpy.array(p1) + trans
    car_p2 = frac_to_car(lattice, p2)
    for each_p1 in all_p1:
        car_p1 = frac_to_car(lattice, each_p1)
        if radius - distance(car_p1, car_p2) > symprec:
            return True
            break
    return False

def _is_in_ele(ele, l_ele):
    if l_ele is None:
        # 若list_ele列表为None则表示删除所有找到的位点
        return True
    else:
        return ele in l_ele
