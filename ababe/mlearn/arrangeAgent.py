# coding: utf-8
# Distributed under the terms of the MIT License.

import ababe.mlearn.arrange as arrange
from ababe.stru.scaffold2d import SitesGrid2d, CStru2d
from ababe.stru.scaffold import SitesGrid, CStru
from ababe.stru.element import GhostSpecie

class Boron2dArrangeProblem(arrange.ArrangeProblem):

    def __init__(self, goal_num):
        self.goal_num = goal_num
        m = [[2, 0], [1, 1.732058]]
        self.delta = 2.1    # the distance threshold to collect neighbors
        sg = SitesGrid2d([[GhostSpecie() for _ in range(goal_num)]
                                            for _ in range(goal_num)])
        self.bg = CStru2d(m, sg)
        self.start_position = self.bg.get_midpoint()

    def get_start_position(self):
        return self.start_position

    def is_goal(self, len_path):
        is_goal = len_path == self.goal_num

        return is_goal

    def get_neighbors(self, position):
        # neighbors is a list contain positions
        neighbors = self.bg.get_neighbors(position, self.delta)

        return neighbors

class FCC3dArrangeProblem(arrange.ArrangeProblem):

    def __init__(self, goal_num):
        self.goal_num = goal_num
        m =         [[0, 0.5, 0.5],
                    [0.5, 0, 0.5],
                    [0.5, 0.5, 0]]
        self.delta = 2.1    # the distance threshold to collect neighbors
        sg = SitesGrid([[[GhostSpecie() for _ in range(goal_num)]
                                            for _ in range(goal_num)]
                                                for _ in range(goal_num)])
        self.bg = CStru(m, sg)
        self.start_position = self.bg.get_midpoint()

    def get_start_position(self):
        return self.start_position

    def is_goal(self, len_path):
        is_goal = len_path == self.goal_num

        return is_goal

    def get_neighbors(self, position):
        # neighbors is a list contain positions
        neighbors = self.bg.get_neighbors(position, self.delta)

        return neighbors
