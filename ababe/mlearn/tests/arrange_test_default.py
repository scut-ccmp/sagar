# coding: utf-8
# Distributed under the terms of the MIT License.

import nose
from nose.tools import *

from ababe.mlearn.arrangeAgent import Boron2dArrangeProblem, FCC3dArrangeProblem
import ababe.mlearn.arrange as arrange

from ababe.stru.scaffold2d import SitesGrid2d, CStru2d
from ababe.stru.element import GhostSpecie

class TestBoronArrange(object):

    def setUp(self):
        self.problom = Boron2dArrangeProblem(50)

    def test_path_result(self):
        path = arrange.random_walk_arrange(self.problom)
        eq_(len(set(path)), 50)

class TestFCCArrange(object):

    def setUp(self):
        self.problom = FCC3dArrangeProblem(20)

    def test_path_result(self):
        path = arrange.random_walk_arrange(self.problom)
        eq_(len(set(path)), 20)
        # eq_(set(path), set[(1,0)])
if __name__ == "__main__":
    nose.main()