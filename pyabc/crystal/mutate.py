# 该模块包含针对`pyabc.crystal.structure`中
# 的MutableCell产生的对象的操作。
import numpy

from pyabc.crystal.structure import MutableCell
from pyabc.crystal.structure import get_symbol, car_to_frac, frac_to_car

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

def remove_sites_in_a_circle(mcell, cc, radius, ele=None):
    """
    移除以cc为圆心，radius为半径的圆内，所有ele的原子。
    cc: 圆心 (x, y, z) 笛卡尔坐标
    radius: 半径 单位为A
    ele: 要删除的元素的列表，若为None，则删除所有元素找到的位点
    """
    pass
