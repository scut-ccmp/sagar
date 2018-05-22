import unittest
import numpy
import copy

from pyabc.crystal.structure import Cell
from pyabc.crystal.utils import non_dup_hnfs, _is_hnf_dup, _hnfs
from pyabc.crystal.utils import IntMat3x3, extended_gcd, snf, atoms_gen
from pyabc.crystal.utils import HartForcadePermutationGroup as HFPG


class TestHnf(unittest.TestCase):

    def test(self):
        pass

    def setUp(self):
        # BCC
        bcc_latt = [1, 1, -1,
                    -1, 1, 1,
                    1, -1, 1]
        bcc_pos = [(0, 0, 0)]
        bcc_atoms = [0]
        self.bcc_pcell = Cell(bcc_latt, bcc_pos, bcc_atoms)

        # FCC
        fcc_latt = [0, 5, 5,
                    5, 0, 5,
                    5, 5, 0]
        fcc_pos = [(0, 0, 0)]
        fcc_atoms = [0]
        self.fcc_pcell = Cell(fcc_latt, fcc_pos, fcc_atoms)

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
        got = [len(non_dup_hnfs(self.bcc_pcell, i))
               for i in range(1, 8)]
        # for h in non_dup_hnfs(self.bcc_pcell, 7):
        #     print(h)
        self.assertEqual(got, wanted)

        # FCC
        wanted = [1, 2, 3, 7, 5, 10, 7]
        got = [len(non_dup_hnfs(self.fcc_pcell, i))
               for i in range(1, 8)]
        self.assertEqual(got, wanted)

        # HCP
        wanted = [1, 3, 5, 11, 7, 19, 11, 34]
        got = [len(non_dup_hnfs(self.hcp_pcell, i))
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

        # debug for compare method
        # numpy.mod problem!
        hnf_x = numpy.array([[1, 0, 0],
                             [0, 1, 2],
                             [0, 0, 5]])
        hnf_y = numpy.array([[1, 0, 3],
                             [0, 1, 3],
                             [0, 0, 5]])
        rot_syms = self.bcc_pcell.get_rotations(1e-3)

        is_dup = _is_hnf_dup(hnf_x, hnf_y, rot_syms, prec=1e-5)
        self.assertTrue(is_dup)

        # debug for compare method
        # numpy.astype problem!
        hnf_x = numpy.array([[1, 0, 6],
                             [0, 1, 6],
                             [0, 0, 7]])
        hnf_y = numpy.array([[1, 0, 3],
                             [0, 1, 6],
                             [0, 0, 7]])
        rot_syms = self.bcc_pcell.get_rotations(1e-3)

        is_dup = _is_hnf_dup(hnf_x, hnf_y, rot_syms, prec=1e-5)
        self.assertTrue(is_dup)


class TestMat3x3(unittest.TestCase):

    def setUp(self):
        self.mat = IntMat3x3([0, 1, 2,
                              3, 4, 5,
                              6, 7, 8])
        self.realmat = IntMat3x3([2, 4, 4,
                                  -6, 6, 12,
                                  10, -4, -16])

    def test_get_snf(self):
        mat = copy.copy(self.realmat)
        snf_D, snf_S, snf_T = mat.get_snf()
        SAT = numpy.matmul(snf_S, numpy.matmul(self.realmat.mat, snf_T))

        wanted_mat = numpy.array([2, 0, 0,
                                  0, 6, 0,
                                  0, 0, 12]).reshape((3, 3))
        numpy.testing.assert_almost_equal(SAT, wanted_mat)
        numpy.testing.assert_almost_equal(snf_D, wanted_mat)

    def test_dead_loop_bug(self):
        mat = IntMat3x3([1, 0, 0,
                         1, 1, 0,
                         0, 0, 7])
        # import pdb; pdb.set_trace()
        ori_mat = copy.copy(mat)
        snf_D, snf_S, snf_T = mat.get_snf()
        SAT = numpy.matmul(snf_S, numpy.matmul(ori_mat.mat, snf_T))

        wanted_mat = numpy.array([1, 0, 0,
                                  0, 1, 0,
                                  0, 0, 7]).reshape((3, 3))
        numpy.testing.assert_almost_equal(SAT, wanted_mat)
        numpy.testing.assert_almost_equal(snf_D, wanted_mat)

    def test_diag_increment_bug(self):
        mat = IntMat3x3([1, 0, 0,
                         0, 2, 0,
                         0, 0, 1])
        ori_mat = copy.copy(mat)
        snf_D, snf_S, snf_T = mat.get_snf()
        SAT = numpy.matmul(snf_S, numpy.matmul(ori_mat.mat, snf_T))

        wanted_mat = numpy.array([1, 0, 0,
                                  0, 1, 0,
                                  0, 0, 2]).reshape((3, 3))
        numpy.testing.assert_almost_equal(SAT, wanted_mat)
        numpy.testing.assert_almost_equal(snf_D, wanted_mat)

    def test_snf_random(self):
        for i in range(100):
            mat = numpy.random.randint(100, size=9).reshape((3, 3))
            mat = IntMat3x3(mat)
            ori_mat = copy.copy(mat)
            snf_D, snf_S, snf_T = mat.get_snf()
            SAT = numpy.matmul(snf_S, numpy.matmul(ori_mat.mat, snf_T))
            numpy.testing.assert_almost_equal(SAT, snf_D)

    def test_snf_diag(self):
        for i in range(100):
            mat = numpy.random.randint(100, size=9).reshape((3, 3))
            mat = IntMat3x3(mat)
            snf_D, _, _ = mat.get_snf()
            self.assertTrue(self.is_diag(snf_D))

    def is_diag(self, mat):
        return numpy.all(mat == numpy.diag(numpy.diagonal(mat)))

    def test_snf_diag_positive(self):
        for i in range(10):
            mat = numpy.random.randint(100, size=9).reshape((3, 3))
            mat = IntMat3x3(mat)
            self.assertTrue(numpy.all(mat.mat >= numpy.zeros_like(mat.mat)))

    def test_snf_diag_incremental(self):
        for i in range(10):
            mat = numpy.random.randint(100, size=9).reshape((3, 3))
            mat = IntMat3x3(mat)
            list_diag = numpy.diagonal(mat.mat).tolist()
            self.assertTrue(sorted(list_diag), list_diag)

    def test_search_first_pivot(self):
        self.assertEqual(self.mat.search_first_pivot(), 1)

    def test_swap_rows(self):
        mat = copy.copy(self.mat)
        mat.swap_rows(0, 1)
        wanted_mat = numpy.array([3, 4, 5,
                                  0, 1, 2,
                                  6, 7, 8]).reshape((3, 3))
        wanted_op = numpy.array([0, 1, 0,
                                 1, 0, 0,
                                 0, 0, 1]).reshape((3, 3))
        numpy.testing.assert_almost_equal(mat.mat, wanted_mat)
        numpy.testing.assert_almost_equal(mat.opL, wanted_op)

        # make sure operation is right, which can restore origin matrix
        wanted_ori_mat = self.mat.mat
        got = numpy.matmul(mat.opL, numpy.matmul(mat.mat, mat.opR))
        numpy.testing.assert_almost_equal(got, wanted_ori_mat)

    def test_flip_sign_row(self):
        mat = copy.copy(self.mat)
        mat.flip_sign_row(1)
        wanted_mat = numpy.array([0, 1, 2,
                                  -3, -4, -5,
                                  6, 7, 8]).reshape((3, 3))
        wanted_op = numpy.array([1, 0, 0,
                                 0, -1, 0,
                                 0, 0, 1]).reshape((3, 3))
        numpy.testing.assert_almost_equal(mat.mat, wanted_mat)
        numpy.testing.assert_almost_equal(mat.opL, wanted_op)

        # make sure operation is right, which can restore origin matrix
        wanted = mat.mat
        got = numpy.matmul(mat.opL, numpy.matmul(self.mat.mat, mat.opR))
        numpy.testing.assert_almost_equal(got, wanted)

    def test_set_zero(self):
        mat = IntMat3x3([3, 4, 5,
                         0, 1, 2,
                         6, 7, 8])
        ori_mat = copy.copy(mat)
        r, s, t = extended_gcd(mat.mat[0, 0], mat.mat[2, 0])
        mat._set_zero(0, 2, mat.mat[0, 0], mat.mat[2, 0], r, s, t)
        wanted_mat = numpy.array([3, 4, 5,
                                  0, 1, 2,
                                  0, -1, -2]).reshape((3, 3))
        wanted_op = numpy.array([1, 0, 0,
                                 0, 1, 0,
                                 -2, 0, 1]).reshape((3, 3))
        numpy.testing.assert_almost_equal(mat.mat, wanted_mat)
        numpy.testing.assert_almost_equal(mat.opL, wanted_op)

        # make sure operation is right, which can restore origin matrix
        wanted = mat.mat
        got = numpy.matmul(mat.opL, numpy.matmul(ori_mat.mat, mat.opR))
        numpy.testing.assert_almost_equal(got, wanted)

    def test_zero_first_column(self):
        mat = copy.copy(self.realmat)
        mat._zero_first_column()
        wanted_mat = numpy.array([2, 4, 4,
                                  0, 18, 24,
                                  0, -24, -36]).reshape((3, 3))
        wanted_op = numpy.array([1, 0, 0,
                                 3, 1, 0,
                                 -5, 0, 1]).reshape((3, 3))
        numpy.testing.assert_almost_equal(mat.mat, wanted_mat)
        numpy.testing.assert_almost_equal(mat.opL, wanted_op)

        # make sure operation is right, which can restore origin matrix
        wanted = mat.mat
        got = numpy.matmul(mat.opL, numpy.matmul(self.realmat.mat, mat.opR))
        numpy.testing.assert_almost_equal(got, wanted)

    def test_zero_first_ele_in_row_i(self):
        mat = copy.copy(self.realmat)
        mat._zero_first_ele_in_row_i(1)
        wanted_mat = numpy.array([2, 4, 4,
                                  0, -18, -24,
                                  10, -4, -16]).reshape((3, 3))
        wanted_op = numpy.array([1, 0, 0,
                                 -3, -1, 0,
                                 0, 0, 1]).reshape((3, 3))
        numpy.testing.assert_almost_equal(mat.mat, wanted_mat)
        numpy.testing.assert_almost_equal(mat.opL, wanted_op)

        # make sure operation is right, which can restore origin matrix
        wanted = mat.mat
        got = numpy.matmul(mat.opL, numpy.matmul(self.realmat.mat, mat.opR))
        numpy.testing.assert_almost_equal(got, wanted)

    def test_first_exact_division(self):
        mat = IntMat3x3([1, 0, 0,
                         1, 1, 0,
                         0, 0, 7])
        ori_mat = copy.copy(mat)
        mat._first_exact_division()
        wanted_mat = numpy.array([1, 0, 0,
                                  0, 1, 0,
                                  0, 0, 7]).reshape((3, 3))
        wanted_op = numpy.array([1, 0, 0,
                                 -1, 1, 0,
                                 0, 0, 1]).reshape((3, 3))
        numpy.testing.assert_almost_equal(mat.mat, wanted_mat)
        numpy.testing.assert_almost_equal(mat.opL, wanted_op)

        # make sure operation is right, which can restore origin matrix
        wanted = mat.mat
        got = numpy.matmul(mat.opL, numpy.matmul(ori_mat.mat, mat.opR))
        numpy.testing.assert_almost_equal(got, wanted)

    def test_zero_first_row(self):
        mat = copy.copy(self.realmat)
        mat._zero_first_row()
        wanted_mat = numpy.array([2, 0, 0,
                                  -6, 18, 24,
                                  10, -24, -36]).reshape((3, 3))
        wanted_op = numpy.array([1, -2, -2,
                                 0, 1, 0,
                                 0, 0, 1]).reshape((3, 3))
        numpy.testing.assert_almost_equal(mat.mat, wanted_mat)
        numpy.testing.assert_almost_equal(mat.opR, wanted_op)

        # make sure operation is right, which can restore origin matrix
        wanted = mat.mat
        got = numpy.matmul(mat.opL, numpy.matmul(self.realmat.mat, mat.opR))
        numpy.testing.assert_almost_equal(got, wanted)

    def test_zero_second_column(self):
        mat = IntMat3x3([2, 0, 0,
                         0, 6, 12,
                         0, 18, 24])
        ori_mat = copy.copy(mat)
        mat._zero_second_column()
        wanted_mat = numpy.array([2, 0, 0,
                                  0, 6, 12,
                                  0, 0, -12]).reshape((3, 3))
        wanted_op = numpy.array([1, 0, 0,
                                 0, 1, 0,
                                 0, -3, 1]).reshape((3, 3))
        numpy.testing.assert_almost_equal(mat.mat, wanted_mat)
        numpy.testing.assert_almost_equal(mat.opL, wanted_op)

        # make sure operation is right, which can restore origin matrix
        wanted = mat.mat
        got = numpy.matmul(mat.opL, numpy.matmul(ori_mat.mat, mat.opR))
        numpy.testing.assert_almost_equal(got, wanted)

    def test_second_exact_division(self):
        mat = IntMat3x3([1, 0, 0,
                         0, 1, 0,
                         0, 1, 7])
        ori_mat = copy.copy(mat)
        mat._second_exact_division()
        wanted_mat = numpy.array([1, 0, 0,
                                  0, 1, 0,
                                  0, 0, 7]).reshape((3, 3))
        wanted_op = numpy.array([1, 0, 0,
                                 0, 1, 0,
                                 0, -1, 1]).reshape((3, 3))
        numpy.testing.assert_almost_equal(mat.mat, wanted_mat)
        numpy.testing.assert_almost_equal(mat.opL, wanted_op)

        # make sure operation is right, which can restore origin matrix
        wanted = mat.mat
        got = numpy.matmul(mat.opL, numpy.matmul(ori_mat.mat, mat.opR))
        numpy.testing.assert_almost_equal(got, wanted)

    def test_zero_second_row(self):
        mat = IntMat3x3([2, 0, 0,
                         0, 6, 12,
                         0, 0, -12])
        ori_mat = copy.copy(mat)
        mat._zero_second_row()
        wanted_mat = numpy.array([2, 0, 0,
                                  0, 6, 0,
                                  0, 0, -12]).reshape((3, 3))
        wanted_op = numpy.array([1, 0, 0,
                                 0, 1, -2,
                                 0, 0, 1]).reshape((3, 3))
        numpy.testing.assert_almost_equal(mat.mat, wanted_mat)
        numpy.testing.assert_almost_equal(mat.opR, wanted_op)

        # make sure operation is right, which can restore origin matrix
        wanted = mat.mat
        got = numpy.matmul(mat.opL, numpy.matmul(ori_mat.mat, mat.opR))
        numpy.testing.assert_almost_equal(got, wanted)

    def test_positive_mat_2_2(self):
        mat = IntMat3x3([2, 0, 0,
                         0, 6, 0,
                         0, 0, -12])
        ori_mat = copy.copy(mat)
        mat._positive_mat_2_2()
        wanted_mat = numpy.array([2, 0, 0,
                                  0, 6, 0,
                                  0, 0, 12]).reshape((3, 3))
        wanted_op = numpy.array([1, 0, 0,
                                 0, 1, 0,
                                 0, 0, -1]).reshape((3, 3))
        numpy.testing.assert_almost_equal(mat.mat, wanted_mat)
        numpy.testing.assert_almost_equal(mat.opL, wanted_op)

        # make sure operation is right, which can restore origin matrix
        wanted = mat.mat
        got = numpy.matmul(mat.opL, numpy.matmul(ori_mat.mat, mat.opR))
        numpy.testing.assert_almost_equal(got, wanted)

    def test_sort_diag(self):
        mat = IntMat3x3([1, 0, 0,
                         0, 2, 0,
                         0, 0, 1])
        ori_mat = copy.copy(mat)
        mat._sort_diag()
        wanted_mat = numpy.array([1, 0, 0,
                                  0, 1, 0,
                                  0, 0, 2]).reshape((3, 3))
        wanted_opL = numpy.array([1, 0, 0,
                                  0, 0, 1,
                                  0, 1, 0]).reshape((3, 3))
        wanted_opR = numpy.array([1, 0, 0,
                                  0, 0, 1,
                                  0, 1, 0]).reshape((3, 3))
        numpy.testing.assert_almost_equal(mat.mat, wanted_mat)
        numpy.testing.assert_almost_equal(mat.opL, wanted_opL)
        numpy.testing.assert_almost_equal(mat.opR, wanted_opR)

        # make sure operation is right, which can restore origin matrix
        wanted = mat.mat
        got = numpy.matmul(mat.opL, numpy.matmul(ori_mat.mat, mat.opR))
        numpy.testing.assert_almost_equal(got, wanted)


class TestSnfHnf(unittest.TestCase):

    def setUp(self):
        # BCC
        bcc_latt = [0.5, 0.5, -0.5,
                    -0.5, 0.5, 0.5,
                    0.5, -0.5, 0.5]
        bcc_pos = [(0, 0, 0)]
        bcc_atoms = [0]
        self.bcc_pcell = Cell(bcc_latt, bcc_pos, bcc_atoms)

        # FCC
        fcc_latt = [0, 5, 5,
                    5, 0, 5,
                    5, 5, 0]
        fcc_pos = [(0, 0, 0)]
        fcc_atoms = [0]
        self.fcc_pcell = Cell(fcc_latt, fcc_pos, fcc_atoms)

    def test_hart_forcade_2008_table_III(self):
        """
        TAKE CARE! The second line of table is snfs of hnfs which are non-redundant
        """
        wanted_a = [1, 7, 13, 35, 31, 91, 57, 155,
                    130, 217, 133, 455, 183, 399, 403, 651]
        # 此为hnf去除旋转对称性后在做snf的结果！
        wanted_b = [1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 2, 1, 1, 1, 4]
        wanted_b_quick = [1, 1, 1, 2, 1, 2, 1, 3, 2, 2, 1, 4, 1, 2, 2, 4]

        # duplicated hnfs produce test
        a = []
        for i in range(1, 17):
            len_volume = len([h for h in _hnfs(i)])
            a.append(len_volume)
        self.assertEqual(a, wanted_a)

        # non-duplicated snfs: b slow test 258s
        b = []
        for i in range(1, 17):
            s_set = set()
            for h in non_dup_hnfs(self.fcc_pcell, volume=i):
                snf_D, _, _ = snf(h)
                s_flat_tuple = tuple(numpy.diagonal(snf_D).tolist())
                s_set.add(s_flat_tuple)
            b.append(len(s_set))
        self.assertEqual(b, wanted_b)

        # duplicated snfs: b quick test
        b = []
        for i in range(1, 17):
            s_set = set()
            for h in _hnfs(i):
                snf_D, _, _ = snf(h)
                s_flat_tuple = tuple(numpy.diagonal(snf_D).tolist())
                s_set.add(s_flat_tuple)
            b.append(len(s_set))
        self.assertEqual(b, wanted_b_quick)


class TestHFPG(object):
    """
    有关该类的测试均为输出型测试，因都是中间件，所有不容易对输出有一个良好的正确的评估。
    在使用时取消各个函数的注释来查看函数的输出，以理解函数的行为和作用。
    """

    def setUp(self):
        # FCC
        fcc_latt = [0, 5, 5,
                    5, 0, 5,
                    5, 5, 0]
        fcc_pos = [(0, 0, 0)]
        fcc_atoms = [0]
        self.fcc_pcell = Cell(fcc_latt, fcc_pos, fcc_atoms)

    def test_get_pure_translations(self):
        # fcc_latt = [0, 5, 5,
        #             5, 0, 5,
        #             5, 5, 0]
        # fcc_pos = [(0, 0, 0),
        #            (0.5,0.5,0.5)]
        # fcc_atoms = [0, 2]
        # fcc_pcell = Cell(fcc_latt, fcc_pos, fcc_atoms)
        # hnfs = non_dup_hnfs(fcc_pcell, 4)
        # h = hnfs[4]
        # print("hnf is:")
        # hfpg = HFPG(fcc_pcell, h)
        # print("pure translations is:")
        # print(hfpg.get_pure_translations())
        # print('\n')
        # pass

        # hnfs = non_dup_hnfs(self.fcc_pcell, 8)
        # for h in hnfs:
        #     print("hnf is:")
        #     print(h)
        #     print("snf is:")
        #     snf_D, _, _ = snf(h)
        #     print(snf_D)
        #     hfpg = HFPG(self.fcc_pcell, h)
        #     print("pure translation is:")
        #     print(hfpg.get_pure_translations())
        #     print('\n')
        pass

    def test_exchanged_new_labels(self):
        # hnfs = non_dup_hnfs(self.fcc_pcell, 4)
        # for h in hnfs:
        #     print("hnf is:")
        #     print(h)
        #     print("super cell atoms is:")
        #     print(self.fcc_pcell.extend(h).atoms)
        #     hfpg = HFPG(self.fcc_pcell, h)
        #     print("label exchange is:")
        #     print(hfpg.get_exchanged_new_labels())
        pass

    def test_get_pure_rotations(self):
        # fcc_latt = [0, 5, 5,
        #             5, 0, 5,
        #             5, 5, 0]
        # fcc_pos = [(0, 0, 0),
        #            (0.5,0.5,0.5)]
        # fcc_atoms = [0, 2]
        # fcc_pcell = Cell(fcc_latt, fcc_pos, fcc_atoms)
        # hnfs = non_dup_hnfs(fcc_pcell, 4)
        # h = hnfs[4]
        # print("hnf is:")
        # print(h)
        # hfpg = HFPG(fcc_pcell, h)
        # print("pure rotation is:")
        # print(hfpg.get_pure_rotations())
        # print('\n')
        pass

    def test_get_symmetry(self):
        # fcc_latt = [0, 5, 5,
        #             5, 0, 5,
        #             5, 5, 0]
        # fcc_pos = [(0, 0, 0),
        #            (0.5,0.5,0.5)]
        # fcc_atoms = [0, 2]
        # fcc_pcell = Cell(fcc_latt, fcc_pos, fcc_atoms)
        # hnfs = non_dup_hnfs(fcc_pcell, 4)
        # h = hnfs[4]
        # print("hnf is:")
        # hfpg = HFPG(fcc_pcell, h)
        # print("perm is:")
        # print(hfpg.get_symmetry())
        # print('\n')
        pass


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

    def test_atoms_gen(self):
        input = [2, 2, 2]
        wanted = [(0, 0, 0),
                  (0, 0, 1),
                  (0, 1, 0),
                  (0, 1, 1),
                  (1, 0, 0),
                  (1, 0, 1),
                  (1, 1, 0),
                  (1, 1, 1)]
        for i, got in enumerate(atoms_gen(input)):
            self.assertEqual(got, wanted[i])


if __name__ == "__main__":
    import nose2
    nose2.main()
