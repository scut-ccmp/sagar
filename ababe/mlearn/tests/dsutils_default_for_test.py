# coding: utf-8
# Distributed under the terms of the MIT License.

import random
import unittest

from ababe.mlearn.dsutils import RandomPopSet

# from unittest.mock import patch

# @patch('random.randint', return_value=3)
class RandomPopSet(unittest.TestCase):

    def setUp(self):
        self.rp = RandomPopSet()
        self.rp_nn = RandomPopSet([2,4,5,5,1])

    def test_add(self):
        rp = RandomPopSet()
        rp.add(5)
        rp.add(5)
        self.assertEqual(rp, RandomPopSet([5]))
        rp.add(34)
        self.assertEqual(rp, RandomPopSet([5, 34]))

    def test_len(self):
        self.assertEqual(len(self.rp_nn), 4)
        self.assertEqual(len(self.rp), 0)

    def test_discard(self):
        rp_nn = RandomPopSet([2,4,5,5,1])
        rp_nn.discard(2)
        self.assertEqual(rp_nn, RandomPopSet([5,5,5,1,4]))

    def test_random_pop(self):
        rp_nn = RandomPopSet([2,4,5,5,1])
        random.seed(0)
        self.assertEqual(rp_nn.random_pop(), 1)
        self.assertEqual(rp_nn, RandomPopSet([2,4,5]))

if __name__ == "__main__":
    import nose2
    nose2.main()
