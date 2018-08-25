# 该模块包含针对`pyabc.crystal.structure`中
# 的MutableCell产生的对象的操作。
import numpy

from pyabc.crystal.structure import MutableCell
from pyabc.crystal.structure import get_symbol, car_to_frac, frac_to_car
from pyabc.utils.math import distance


def cell_to_mcell(cell):
    lattice = cell.lattice
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
    cc: 圆心 (x, y, z) 笛卡尔坐标
    radius: 半径 单位为A
    ele: 要删除的元素的list，若为None，则删除所有元素找到的位点
    """
    # TODO: 不用copy
    import copy
    for s in copy.deepcopy(mcell._sites):
        car_cc = frac_to_car(mcell._lattice, cc)
        car_site = frac_to_car(mcell._lattice, s[0])
        if _is_close_in_radius(car_cc, car_site, radius) and _is_in_ele(s[1], list_ele):
            mcell._sites.remove(s)

def _is_close_in_radius(p1, p2, radius):
    # 截断圆对应的半径大于距离，说明两点在圆内，需要删除：返回True
    symprec = 1e-5
    return radius - distance(p1, p2) > symprec

def _is_in_ele(ele, l_ele):
    if l_ele is None:
        # 若list_ele列表为None则表示删除所有找到的位点
        return True
    else:
        return ele in l_ele
