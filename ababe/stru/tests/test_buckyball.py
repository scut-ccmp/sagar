# encoding: utf-8
# Distributed under the terms of the MIT License.

import unittest

from ababe.stru.buckyball import Structure
from ababe.stru.element import Specie

class testBuckyballStructure(unittest.TestCase):

    def setUp(self):
        self.bucky = Structure([6]*60)
        self.speckle = Specie("B")

    def test_a(self):
        pass

    def test_all_speckle_gen(self):
        # all_nondup_structure = Structure.all_speckle_gen(self.bucky, 3, self.speckle)
        # self.assertIn(a, all_nondup_structure)
        init_buckyball = Structure([6]*60)
        gen_with_3speckle = init_buckyball.all_speckle_gen(init_buckyball, 3, self.speckle)

        l = []
        for nodup_gen in gen_with_3speckle:
            l.append(len([i for i in nodup_gen]))
        self.assertEqual(l, [1,1,23,303])

    def test_add_one_speckle(self):
        gen = self.bucky.to_gen()
        pass

    def test_time_of_speckle_gen(self):
        import time
        init_buckyball = Structure([6]*60)
        gen_with_5speckle = init_buckyball.all_speckle_gen(init_buckyball, 3, self.speckle)

        for nodup_gen in gen_with_5speckle:
            start = time.time()
            _ = [i for i in nodup_gen]
            end = time.time()
            t_spend = end - start
            print(t_spend)

        self.assertLess(t_spend, 1)

if __name__ == "__main__":
    import nose2
    nose2.main()
