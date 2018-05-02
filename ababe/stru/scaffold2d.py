# coding: utf-8
# Distributed under the terms of the MIT License.
from ababe.stru.element import Specie, GhostSpecie
import numpy as np

from scipy.spatial import cKDTree
from operator import itemgetter

class SitesGrid2d(object):

    def __init__(self, sites):
        self._sites = sites
        self._width = len(sites)
        self._length = len(sites[0])

    @property
    def width(self):
        return self._width

    @property
    def length(self):
        return self._length

    @classmethod
    def sea(cls, width, length, sp = GhostSpecie()):
        sites = [[sp for _ in range(length)] 
                        for _ in range(width)]

        return cls(sites)

    def __eq__(self, other):
        if other == None: return False
        return self._sites == other._sites

    def get_array(self):
        mfunc = np.vectorize(lambda sp: sp.Z)
        arr = mfunc(np.array(self._sites))
        return arr

    @classmethod
    def from_array(cls, arr):
        mfunc = np.vectorize(lambda n: Specie.to_sp(n))
        sarr = mfunc(arr)
        return cls(sarr.tolist())

class CStru2d(object):
    def __init__(self, m, sg):
        self._matrix = m
        self.width = sg.width
        self.length = sg.length
        self._sites_grid = sg

    @property
    def m(self):
        return self._matrix

    def get_array(self):
        return self._sites_grid.get_array()

    def get_midpoint(self):
        w = self.width
        l = self.length
        return (w//2, l//2)

    # @staticmethod
    # def _pos2coor(pos):
    #     a, b = np.array(self.m)
    #     x, y = pos
    #     coor = a*x + b*y    # an array
    #     return tuple(coor)

    def get_neighbors(self, pos, delta):
        
        def _pos2coor(pos):
            a, b = np.array(self.m)
            x, y = pos
            coor = a*x + b*y    # an array
            return tuple(coor)

        def p_gen():
            for x in range(self.width):
                for y in range(self.length):
                    yield(x, y)
        
        point = _pos2coor(pos)
        # w = self.width
        # l = self.length
        coor_map = {p : _pos2coor(p) for p in p_gen()}
        del coor_map[pos]

        points = list(coor_map.values())
        points_tree = cKDTree(points)
        ind = points_tree.query_ball_point(point, delta)
        neighbors = itemgetter(*ind)(list(coor_map.keys()))
        return set(neighbors)