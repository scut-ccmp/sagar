import unittest
import numpy

from pyyabc.utils.math import distance, closest_pair, extended_gcd


class TestCommonUtils(unittest.TestCase):

    def test_extended_gcd(self):
        self._test_extended_gcd_n_times(100)

    def _test_extended_gcd_n_times(self, n):
        for i in range(n):
            aa, bb = numpy.random.randint(10, size=2) + 1
            r, s, t = extended_gcd(aa, bb)
            # print("%d = %d * (%d) + %d * (%d)" %
            #       (r, aa, s, bb, t))
            wanted = aa * s + bb * t
            self.assertEqual(r, wanted)

    def test_is_int_np_array(self):
        pass

    def test_binomial_coeff(self):
        pass

    def test_refine_positions(self):
        pass

    def test_distance(self):
        p1 = (0, 0, 0)
        p2 = (0, 3, 4)
        self.assertEqual(distance(p1, p2), 5.0)

    def test_closest_pair(self):
        points = numpy.array([0, 0, 0,
                              0, 3, 4,
                              3, 60, 6,
                              10.1, 10.2, 10.3]).reshape((-1, 3))
        self.assertEqual(closest_pair(points), 5)


if __name__ == "__main__":
    import nose2
    nose2.main()
