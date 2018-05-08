import unittest
import numpy

from pyabc.crystal.structure import Cell
from pyabc.crystal.derive import hnf_cells, _is_hnf_dup
from pyabc.crystal.derive import snf, extended_gcd

from pyabc.crystal.derive import _search_first_pivot, _swap_rows, \
    _flip_sign_row, _set_zero, _zero_first_column


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

    def setUp(self):
        self.mat = numpy.array([0, 1, 2,
                                3, 4, 5,
                                6, 7, 8]).reshape((3, 3))
        self.realmat = numpy.array([2, 4, 4,
                                    -6, 6, 12,
                                    10, -4, -16]).reshape((3, 3))

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
        mat = self.realmat
        snf_S, snf_A, snf_T = snf(mat)
        SAT = numpy.matmul(snf_S, numpy.matmul(mat, snf_T))

        wanted_mat = numpy.array([2, 0, 0,
                                  0, 6, 0,
                                  0, 0, 12]).reshape((3, 3))

        self.assertTrue(numpy.allclose(SAT, wanted_mat))
        self.assertTrue(numpy.allclose(snf_A, wanted_mat))

    # strange but robust test
    def test_smith_normal_form_random(self):
        self._test_Instantiate(100)

    def _test_Instantiate(self, n):
        for i in range(n):
            mat = self._get_random_mat()
            snf_S, snf_A, snf_T = snf(mat)
            SAT = numpy.matmul(snf_S, numpy.matmul(mat, snf_T))

            numpy.testing.assert_almost_equal(numpy.linalg.det(snf_S), 1)
            numpy.testing.assert_almost_equal(numpy.linalg.det(snf_T), 1)

    def _get_random_mat(self):
        k = 15
        mat = numpy.random.randint(k, size=(3, 3)) - k // 2
        if numpy.linalg.det(mat) < 0.5:
            mat = self._get_random_mat()
        return mat

    def test_search_first_pivot(self):
        self.assertEqual(_search_first_pivot(self.mat), 1)

    def test_swap_rows(self):
        mat, op = _swap_rows(self.mat, 0, 1)
        wanted_mat = numpy.array([3, 4, 5,
                                  0, 1, 2,
                                  6, 7, 8]).reshape((3, 3))
        wanted_op = numpy.array([0, 1, 0,
                                 1, 0, 0,
                                 0, 0, 1]).reshape((3, 3))
        numpy.testing.assert_almost_equal(mat, wanted_mat)
        numpy.testing.assert_almost_equal(op, wanted_op)

    def test_flip_sign_row(self):
        mat, op = _flip_sign_row(self.mat, 1)
        wanted_mat = numpy.array([0, 1, 2,
                                  -3, -4, -5,
                                  6, 7, 8]).reshape((3, 3))
        wanted_op = numpy.array([1, 0, 0,
                                 0, -1, 0,
                                 0, 0, 1]).reshape((3, 3))
        numpy.testing.assert_almost_equal(mat, wanted_mat)
        numpy.testing.assert_almost_equal(op, wanted_op)

    def test_set_zero(self):
        mat, _ = _swap_rows(self.mat, 0, 1)
        # now_mat = numpy.array([3, 4, 5,
        #                        0, 1, 2,
        #                        6, 7, 8]).reshape((3, 3))
        r, s, t = extended_gcd(mat[0, 0], mat[2, 0])
        mat, op = _set_zero(mat, 0, 2, mat[0, 0], mat[2, 0], r, s, t)
        wanted_mat = numpy.array([3, 4, 5,
                                  0, 1, 2,
                                  0, -1, -2]).reshape((3, 3))
        wanted_op = numpy.array([1, 0, 0,
                                 0, 1, 0,
                                 -2, 0, 1]).reshape((3, 3))
        numpy.testing.assert_almost_equal(mat, wanted_mat)
        numpy.testing.assert_almost_equal(op, wanted_op)

    def test_zero_first_column(self):
        mat, op = _zero_first_column(self.realmat, 1)
        wanted_mat = numpy.array([2, 4, 4,
                                  0, -18, -24,
                                  10, -4, -16]).reshape((3, 3))
        wanted_op = numpy.array([1, 0, 0,
                                 -3, -1, 0,
                                 0, 0, 1]).reshape((3, 3))
        numpy.testing.assert_almost_equal(mat, wanted_mat)
        numpy.testing.assert_almost_equal(op, wanted_op)


if __name__ == "__main__":
    import nose2
    nose2.main()
