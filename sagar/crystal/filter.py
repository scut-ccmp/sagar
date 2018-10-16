# -*- coding: utf-8 -*-
import numpy as np

from sagar.crystal.structure import Cell
from sagar.toolkit.mathtool import distance


class Filter(object):
    """
    结构过滤器基类
    当返回True时表示该结构满足条件，接受
    当返回False时表示结构不满足条件，不接受
    """
    def is_accepted(self, cell):
        return isinstance(cell, Cell)

class MinDistanceFilter(Filter):
    """
    参数中的element是表示元素的字符串
    最小距离过滤器（连通性过滤器）
    当返回True时表示该结构中指定元素的距离在晶体中不相连，分散，可接受
    当返回False时藐视该结构中指定元素的距离在晶体中过近，相连，不接受
    """
    def __init__(self, element, distance):
        self.element = element
        self.distance = distance

    def is_accepted(self, cell):
        return cell.check(elements = [self.element], limit=self.distance)


class SpaceGroupFilter(Filter):
    """
    对称性过滤器
    当返回True时表示该结构的空间群满足
    当返回False时表示该结构空间群不满足
    """
    def __init__(self, list_spg, symprec=1e-3):
        self.list_spg = [str(x) for x in list_spg]
        self.symprec = symprec

    def is_accepted(self, cell):
        symmetry = cell.get_spacegroup(self.symprec)
        spg_name, spg_num = symmetry.split(' ')
        spg_num = spg_num[1:-1]
        if spg_name in self.list_spg or spg_num in self.list_spg:
            return True
        return False
