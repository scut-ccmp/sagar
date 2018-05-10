import numpy
import copy

from math import sqrt
from itertools import product

from pyabc.crystal.structure import Cell, is_primitive_cell


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


def non_dup_hnfs(pcell, volume=1, symprec=1e-5, comprec=1e-5):
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
                        "please provide pyabc.crystal.structure.Cell object.".format(type(pcell)))

    if not is_primitive_cell(pcell):
        raise ValueError("cell object you provide is not a primitive cell "
                         "Therefore meaningless to get non duplicated hnf cells "
                         "You can use pcell.get_primitive() first.")

    nodup_hnfs = []

    # Using rot without inversion class double speed
    # rot_list = pcell.get_rotations(symprec)
    rot_list = pcell.get_rotations_without_inversion(symprec)
    for hnf in _hnfs(volume):
        if _not_contain(nodup_hnfs, hnf, rot_list, comprec):
            nodup_hnfs.append(hnf)

    return nodup_hnfs


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
        # is_int1 = numpy.all(numpy.isclose(m, m.astype(numpy.int), atol=1e-3))
        # is_int2 = numpy.allclose(numpy.mod(m, 1), numpy.zeros_like(m), atol=prec)
        # is_array_int = numpy.all(numpy.isclose(m, numpy.around(m), atol=prec))
        if is_int_np_array(m, prec):
            return True

    return False


def is_int_np_array(npa, atol=1e-5):
    return numpy.all(numpy.isclose(npa, numpy.around(npa), atol=atol))


def extended_gcd(aa, bb):
    """
    Algorithm: https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Iterative_method_2

    parameters:
    aa, bb: int

    return: r, s, t

    r = s * aa + t * bb
    """
    lastremainder, remainder = abs(aa), abs(bb)
    x, lastx, y, lasty = 0, 1, 1, 0
    while remainder:
        lastremainder, (quotient, remainder) = remainder, divmod(
            lastremainder, remainder)
        x, lastx = lastx - quotient * x, x
        y, lasty = lasty - quotient * y, y
    return lastremainder, lastx * (-1 if aa < 0 else 1), lasty * (-1 if bb < 0 else 1)


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

            self._positive_mat_2_2()

        if not self._is_incremental_diag():
            self._sort_diag()

        return self._mat, self._opL, self._opR

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
        matT = __class__(self._mat.T)
        matT._zero_first_column()
        op = matT.opL.T
        self._mat = numpy.matmul(self._mat, op)
        self._opR = numpy.matmul(self._opR, op)

    def _zero_second_column(self):
        if self._mat[2, 1] % self._mat[1, 1] == 0:
            self._second_exact_division()
        if self._mat[2, 1] < 0:
            self.flip_sign_row(2)
        r, s, t = extended_gcd(self._mat[1, 1], self._mat[2, 1])
        self._set_zero(1, 2, self._mat[1, 1], self._mat[2, 1], r, s, t)

    def _second_exact_division(self):
        op = numpy.eye(3, dtype='int')
        op[2, 1] = -self._mat[2, 1] // self._mat[1, 1]
        self._mat = numpy.matmul(op, self._mat)
        self._opL = numpy.matmul(op, self._opL)

    def _zero_second_row(self):
        matT = __class__(self._mat.T)
        matT._zero_second_column()
        op = matT.opL.T
        self._mat = numpy.matmul(self._mat, op)
        self._opR = numpy.matmul(self._opR, op)

    def _positive_mat_2_2(self):
        if self._mat[2, 2] < 0:
            self.flip_sign_row(2)

    def _sort_diag(self):
        mat_flat = numpy.diagonal(self._mat)
        idx = numpy.argsort(mat_flat)
        op = numpy.eye(3, dtype='int')[idx]
        self._mat = numpy.matmul(op, self._mat)
        self._opL = numpy.matmul(op, self._opL)

        self._mat = numpy.matmul(self._mat, op.T)
        self._opR = numpy.matmul(self._opR, op.T)

# TODO: utils for Hart-Forcade algorithm


class HartForcadePermutationGroup(object):

    def __init__(self, pcell, hnf):
        if not isinstance(pcell, Cell):
            raise TypeError(
                "want pyabc.crystal.structure.Cell, got {:}".format(type(cell)))
        self._cell = pcell
        self._hnf = hnf
        self._snf, _, _ = snf(hnf)
        self._quotient = numpy.diagonal(self._snf).tolist()
        self._volume = numpy.diagonal(self._snf).prod()
        self._nsites = len(pcell.atoms)  # 最小原胞中原子个数 如：hcp为2

    def get_symmetry(self, symprec=1e-5):
        pass

    def get_pure_translations(self, symprec=1e-5):
        itertrans = [list(range(self._quotient[0])),
                     list(range(self._quotient[1])),
                     list(range(self._quotient[2]))]
        size = self._volume
        result = numpy.zeros((size - 1, self._nsites * size), dtype='int') - 1
        iterable = product(*itertrans)

        # remove (0,0,0) 因为它对应的是单位操作，保持原来的构型
        a = next(iterable)  # avoid null translation
        assert a == (0, 0, 0)  # check that first element (0,0,0) is remove

        for t, (i, j, k) in enumerate(iterable):
            for l, m, n in product(*itertrans):
                u = (i + l) % self._quotient[0]
                v = (j + m) % self._quotient[1]
                w = (k + n) % self._quotient[2]
                for s in range(self._nsites):
                    result[t, self._flatten_indices(
                        l, m, n, s)] = self._flatten_indices(u, v, w, s)

        return result

    def _flatten_indices(self, i, j, k, site=0):
        iq = i + site * self._quotient[0]
        jq = j + iq * self._quotient[1]
        kq = k + jq * self._quotient[2]
        return kq

    def get_pure_rotations(self, symprec=1e-5):
        pass

    def get_pure_rotations_without_inversion(self, symprec=1e-5):
        pass
