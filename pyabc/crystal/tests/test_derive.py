import unittest
import numpy

from pyabc.crystal.structure import Cell
from pyabc.crystal.derive import hnf_cells, _is_hnf_dup
from pyabc.crystal.derive import snf, extended_gcd


class TestHnf(unittest.TestCase):

    def test(self):
        pass

    def setUp(self):
        # BCC
        bcc_latt = [0.5, 0.5, -0.5,
                    -0.5, 0.5, 0.5,
                    0.5, -0.5, 0.5]
        bcc_pos = [(0, 0, 0)]
        bcc_atoms = [0]
        self.bcc_pcell = Cell(bcc_latt, bcc_pos, bcc_atoms)

        # HCP
        hcp_b = [2.51900005,  0.,  0.,
                 -1.25950003,  2.18151804, 0.,
                 0., 0.,  4.09100008]
        hcp_positions = [(0.33333334,  0.66666669,  0.25),
                         (0.66666663,  0.33333331,  0.75)]
        hcp_numbers = [0, 0]
        self.hcp_pcell = Cell(hcp_b, hcp_positions, hcp_numbers)

    def test_hnf_cells(self):
        # Results from <PHYSICAL REVIEW B 80, 014120 (2009)>

        # BCC
        wanted = [1, 2, 3, 7, 5, 10, 7]
        got = [len(hnf_cells(self.bcc_pcell, i))
               for i in range(1, 8)]
        self.assertEqual(got, wanted)

        # HCP
        wanted = [1, 3, 5, 11, 7, 19, 11, 34]
        got = [len(hnf_cells(self.hcp_pcell, i))
               for i in range(1, 9)]
        self.assertEqual(got, wanted)

    def test_is_hnf_dup(self):
        hnf_x = numpy.array([[1, 0, 0],
                             [0, 1, 0],
                             [0, 0, 2]])
        hnf_y = numpy.array([[1, 0, 0],
                             [0, 2, 0],
                             [0, 0, 1]])
        rot_syms = self.bcc_pcell.get_rotations(1e-3)
        is_dup = _is_hnf_dup(hnf_x, hnf_y, rot_syms, prec=1e-3)
        self.assertTrue(is_dup)


class TestSnf(unittest.TestCase):

    def test(self):
        pass

    def test_extended_gcd(self):
        self._test_extended_gcd_n_times(100)

    def _test_extended_gcd_n_times(self, n):
        for i in range(n):
            aa, bb = numpy.random.randint(100, size=2) + 1
            r, s, t = extended_gcd(aa, bb)
            # print("%d = %d * (%d) + %d * (%d)" %
            #       (r, aa, s, bb, t))
            wanted = aa * s + bb * t
            self.assertEqual(r, wanted)

    def test_smith_normal_form(self):
        mat = numpy.array([2, 4, 4,
                           -6, 6, 12,
                           10, -4, -16]).reshape((3, 3))
        snf_S, snf_A, snf_T = snf(mat)
        SAT = numpy.dot(snf_S, numpy.dot(mat, snf_T))

        wanted_mat = numpy.array([2, 0, 0,
                                  0, 6, 0,
                                  0, 0, 12]).reshape((3, 3))

        self.assertTrue(numpy.allclose(SAT, wanted_mat))
        self.assertTrue(numpy.allclose(snf_A, wanted_mat))


if __name__ == "__main__":
    import nose2
    nose2.main()
