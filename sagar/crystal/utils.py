# -*- coding: utf-8 -*-
import numpy
import copy

from math import sqrt
from itertools import product

from sagar.crystal.structure import Cell
from sagar.toolkit.mathtool import is_int_np_array, extended_gcd


def _factor(n):
    """
    https://rosettacode.org/wiki/Factors_of_an_integer#Python

    realize that factors come in pairs, the smaller of which is
    no bigger than sqrt(n)
    """
    factors = set()
    for x in range(1, int(sqrt(n)) + 1):
        if n % x == 0:
            factors.add(x)
            factors.add(n // x)
    return sorted(factors)


def _hnfs(det):
    h = []
    for a in _factor(det):
        for d in _factor(det // a):
            f = det // a // d

            for b, c, e in product(range(d), range(f), range(f)):
                yield numpy.array([[a, b, c],
                                   [0, d, e],
                                   [0, 0, f]])


def _hnfs_2D(det):
    # 二维的hnf矩阵生成,z方向为1
    h = []
    for a in _factor(det):
        d = det // a
        for b in range(d):
            yield numpy.array([[a, b, 0],
                               [0, d, 0],
                               [0, 0, 1]])


def non_dup_hnfs(pcell, volume=1, dimension=3, symprec=1e-5, comprec=1e-5):
    """
    hnf_cells return all non duplicated hnf extend cells.

    parameters:

    pcell: Cell object, The primitive cell to be extended
    volume: int, Extend to how large supercellself, default=1
    symprec: int, symmetry precision
    When finding duplicated hnfs the precesion, default=5
    comprec: float, compare precision
    When finding the rotations symmetry of primitive cell, defalut=1e-5

    return:

    A list of 2D numpy.ndarray.
    """
    if not isinstance(pcell, Cell):
        raise TypeError("Can't make hnf cells of {:} "
                        "please provide sagar.crystal.structure.Cell object.".format(type(pcell)))

    if not pcell.is_primitive(symprec):
        raise ValueError("cell object you provide is not a primitive cell "
                         "Therefore meaningless to get non duplicated hnf cells "
                         "You can use pcell.get_primitive() first.")

    nodup_hnfs = []

    # Using rot without inversion class double speed
    # rot_list = pcell.get_rotations(symprec)
    rot_list = pcell.get_rotations_without_inversion(symprec)
    if dimension == 3:
        for hnf in _hnfs(volume):
            if _not_contain(nodup_hnfs, hnf, rot_list, comprec):
                nodup_hnfs.append(hnf)

        return nodup_hnfs
    elif dimension == 2:
        for hnf in _hnfs_2D(volume):
            if _not_contain(nodup_hnfs, hnf, rot_list, comprec):
                nodup_hnfs.append(hnf)

        return nodup_hnfs
    else:
        raise ValueError("Dimension you provide is not as expected."
                         "You can only use dimension = 3(default) or 2.")


def _not_contain(hnf_list, hnf, rot_list, prec):
    for h in hnf_list:
        if _is_hnf_dup(hnf, h, rot_list, prec):
            return False
    return True


def _is_hnf_dup(hnf_x, hnf_y, rot_list, prec=1e-5):
    """
    A hnf act in a cell,
    if H_x * R.T^-1 * H_y ^-1 is an int matrix:
        H_x and H_y produce the same supercell.

    Algorithm: Hermite normal form (HNF) matrices remove translation symmetry duplications
               However, rotation symmetry duplications also need to be removedself.
               That is what this function ``_is_hnf_dup`` do.
    """
    for rot in rot_list:
        m = numpy.matmul(
            numpy.matmul(hnf_x, numpy.linalg.inv(rot.T)),
            numpy.linalg.inv(hnf_y))
        # TODO: stackoverflow add answer
        # is_int1 = numpy.all(numpy.isclose(m, m.astype(numpy.int), atol=1e-3))
        # is_int2 = numpy.allclose(numpy.mod(m, 1), numpy.zeros_like(m), atol=prec)
        # is_array_int = numpy.all(numpy.isclose(m, numpy.around(m), atol=prec))
        if is_int_np_array(m, prec):
            return True

    return False

def snf(mat):
    """
    snf return Smith Normal Form of a 3x3 int matrix
        which result_mat = opL * ori_mat * opR

    parameter:
    mat is a 3x3 or nine int list

    return:
    mat: numpy.ndarray, Smith Normal Form matrix
    opL: left operation
    opR: right operation
    """
    mat = numpy.array(mat).reshape((3, 3))
    intmat = IntMat3x3(mat)
    snf_D, snf_S, snf_T = intmat.get_snf()
    return snf_D, snf_S, snf_T


class IntMat3x3(object):

    def __init__(self, mat):
        # TODO: check if mat not a int matrix
        self._mat = numpy.array(mat, dtype='int').reshape((3, 3))
        self._opL = numpy.eye(3, dtype='int')
        self._opR = numpy.eye(3, dtype='int')

    @property
    def mat(self):
        return self._mat

    @property
    def opL(self):
        return self._opL

    @property
    def opR(self):
        return self._opR

    def is_diag(self):
        return numpy.all(self._mat == numpy.diag(numpy.diagonal(self._mat)))

    def _is_incremental_diag(self):
        is_diag = self.is_diag()
        list_mat_flat = numpy.diagonal(self._mat).tolist()
        is_incremental = sorted(list_mat_flat) == list_mat_flat
        return is_diag and is_incremental

    def search_first_pivot(self):
        for i in range(3):
            if self._mat[i, 0] != 0:
                return i

    def swap_rows(self, i, j):
        """
        return:
        mat: new row permuted mat
        op: is a row permute matrix

        example:

        """
        op = numpy.eye(3, dtype='int')
        op[i, i], op[j, j] = 0, 0
        op[i, j], op[j, i] = 1, 1
        self._mat = numpy.matmul(op, self._mat)
        self._opL = numpy.matmul(op, self._opL)

    def flip_sign_row(self, i):
        """
        mul -1 for i row
        """
        op = numpy.eye(3, dtype='int')
        op[i, i] = -1
        self._mat = numpy.matmul(op, self._mat)
        self._opL = numpy.matmul(op, self._opL)

    def get_snf(self):
        while not self.is_diag():
            pivot = self.search_first_pivot()
            if pivot > 0:
                self.swap_rows(0, pivot)

            self._zero_first_column()
            self._zero_first_row()

            # diagonal lower 2x2 matrix
            self._zero_second_column()
            self._zero_second_row()

            # self._positive_mat_2_2()
            self._positive_diag()

        if not self._is_incremental_diag():
            self._sort_diag()

        if numpy.linalg.det(self._opL) < 0:
            self._opL = self._hand_flip(self._opL)
            self._opR = self._hand_flip(self._opR)

        return self._mat, self._opL, self._opR

    def _hand_flip(self, mat):
        neg = numpy.array([-1, 0, 0,
                           0, -1, 0,
                           0, 0, -1]).reshape((3, 3))
        return numpy.matmul(mat, neg)

    def _set_zero(self, i, j, aa, bb, r, s, t):
        """
        Based on Bezout's identity

        Let mat[i, j] be zero

        return:
        op: is a numpy.ndarray, represent the operations

        example:

        """
        op = numpy.eye(3, dtype='int')
        op[i, i] = s
        op[i, j] = t
        op[j, i] = -bb // r
        op[j, j] = aa // r
        self._mat = numpy.matmul(op, self._mat)
        self._opL = numpy.matmul(op, self._opL)

    def _zero_first_ele_in_row_i(self, i):
        if self._mat[i, 0] < 0:
            self.flip_sign_row(i)
        r, s, t = extended_gcd(self._mat[0, 0], self._mat[i, 0])
        self._set_zero(0, i, self._mat[0, 0], self._mat[i, 0], r, s, t)

    def _zero_first_column(self):
        if self._mat[1, 0] % self._mat[0, 0] == 0 and self._mat[2, 0] % self._mat[0, 0] == 0:
            self._first_exact_division()
        for i in [1, 2]:
            if self._mat[i, 0] != 0:
                self._zero_first_ele_in_row_i(i)

    def _first_exact_division(self):
        op = numpy.eye(3, dtype='int')
        op[1, 0] = -self._mat[1, 0] // self._mat[0, 0]
        op[2, 0] = -self._mat[2, 0] // self._mat[0, 0]
        self._mat = numpy.matmul(op, self._mat)
        self._opL = numpy.matmul(op, self._opL)

    def _zero_first_row(self):
        matT = self.__class__(self._mat.T)
        matT._zero_first_column()
        op = matT.opL.T
        self._mat = numpy.matmul(self._mat, op)
        self._opR = numpy.matmul(self._opR, op)

    def _zero_second_column(self):
        if self._mat[1, 1] == 0 or self._mat[2, 1] % self._mat[1, 1] == 0:
            self._second_exact_division()
        if self._mat[2, 1] < 0:
            self.flip_sign_row(2)
        r, s, t = extended_gcd(self._mat[1, 1], self._mat[2, 1])
        self._set_zero(1, 2, self._mat[1, 1], self._mat[2, 1], r, s, t)

    def _second_exact_division(self):
        op = numpy.eye(3, dtype='int')
        if self._mat[1, 1] == 0:
            op[2, 1] = 0
        else:
            op[2, 1] = -self._mat[2, 1] // self._mat[1, 1]
        self._mat = numpy.matmul(op, self._mat)
        self._opL = numpy.matmul(op, self._opL)

    def _zero_second_row(self):
        matT = self.__class__(self._mat.T)
        matT._zero_second_column()
        op = matT.opL.T
        self._mat = numpy.matmul(self._mat, op)
        self._opR = numpy.matmul(self._opR, op)

    def _positive_diag(self):
        for i in range(0, 3):
            if self._mat[i, i] < 0:
                self.flip_sign_row(i)

    def _sort_diag(self):
        mat_flat = numpy.diagonal(self._mat)
        idx = numpy.argsort(mat_flat)
        op = numpy.eye(3, dtype='int')[idx]
        self._mat = numpy.matmul(op, self._mat)
        self._opL = numpy.matmul(op, self._opL)

        self._mat = numpy.matmul(self._mat, op.T)
        self._opR = numpy.matmul(self._opR, op.T)
