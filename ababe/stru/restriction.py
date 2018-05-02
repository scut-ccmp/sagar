# coding: utf-8
# Distributed under the terms of the MIT License.
from scipy.spatial.distance import pdist
import numpy as np
import re

class Restriction(object):
    """
    The model class for writing specific
    Restriction.
    """

    def __init__(self):
        pass

    def is_satisfied(self, gcell):
        pass

class MinDistanceRestriction(Restriction):
    """
    The class giving the instance of distance constriant.
    Input element specie, distance constriant
    Implet to gcell to check the satisfiction
    Output is whether the structure satisfied the constriant.
    if all distance of ele is > input distance then output True
    else output False
    """

    def __init__(self, tr):
        self.target_ele = tr[0]
        self.target_dist = tr[1]

    def is_satisfied(self, gcell):
        scale = np.array([[2, 0, 0],
                          [0, 2, 0],
                          [0, 0, 2]])
        super_gcell = gcell.supercell(scale)
        target_cart = super_gcell.get_cartesian(ele=self.target_ele)
        # target_cart is a np array of target element's
        # cartesian coordinates
        mindist = np.min(pdist(target_cart))
        is_ok = mindist > self.target_dist
        # import pdb
        # pdb.set_trace()
        return is_ok

class SpaceGroupRestriction(Restriction):
    """
    The class giving the instance of Space group constriant.
    Only input spacegroups get True
    The input is a list of spacegroup No.
    """

    def __init__(self, l_sp):
        self.l_sp = l_sp

    def is_satisfied(self, gcell):
        spg = gcell.get_spacegroup()    # output: 'Pm (6)'
        spg_num = int(re.findall('\d+', spg)[0])
        is_ok = spg_num in self.l_sp

        return is_ok
