import numpy
import numpy as np

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


def hnf_cells(pcell, volume=1, symprec=1e-5, comprec=1e-5):
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

    A list of Cell objects.
    """
    if not isinstance(pcell, Cell):
        raise TypeError("Can't make hnf cells of {:} "
                        "please provide pyabc.crystal.structure.Cell object.".format(type(pcell)))

    if not is_primitive_cell(pcell):
        raise ValueError("cell object you provide is not a primitive cell "
                         "Therefore meaningless to get non duplicated hnf cells "
                         "You can use pcell.get_primitive() first.")

    nodup_hnfs = []
    rot_list = pcell.get_rotations(symprec)
    for hnf in _hnfs(volume):
        if _not_contain(nodup_hnfs, hnf, rot_list, comprec):
            nodup_hnfs.append(hnf)

    nodup_cells = [pcell.extend(hnf) for hnf in nodup_hnfs]
    return nodup_cells


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
    # TODO: efficiency!
    for rot in rot_list:
        m = numpy.matmul(
            numpy.matmul(hnf_x, numpy.linalg.inv(rot.T)),
            numpy.linalg.inv(hnf_y))
        if numpy.allclose(numpy.mod(m, 1), numpy.zeros_like(m), atol=prec):
            return True

    return False


class Snf(object):
    """
    snf 3x3
    """
    def __init__(self, A):
        self._A_orig = np.array(A, dtype='intc')
        self._A = np.array(A, dtype='intc')
        self._Ps = []
        self._Qs = []
        self._L = []
        self._P = None
        self._Q = None
        self._attempt = 0

    def run(self):
        for i in self:
            pass

    def __iter__(self):
        return self

    def __next__(self):
        self._attempt += 1
        if self._first():
            if self._second():
                self._set_PQ()
                raise StopIteration
        return self._attempt

    def next(self):
        self.__next__()

    @property
    def A(self):
        return self._A

    @property
    def P(self):
        return self._P

    @property
    def Q(self):
        return self._Q

    def _set_PQ(self):
        if np.linalg.det(self._A) < 0:
            for i in range(3):
                if self._A[i, i] < 0:
                    self._flip_sign_row(i)
            self._Ps += self._L
            self._L = []

        P = np.eye(3, dtype='intc')
        for _P in self._Ps:
            P = np.dot(_P, P)
        Q = np.eye(3, dtype='intc')
        for _Q in self._Qs:
            Q = np.dot(Q, _Q.T)

        if np.linalg.det(P) < 0:
            P = -P
            Q = -Q

        self._P = P
        self._Q = Q

    def _first(self):
        self._first_one_loop()
        A = self._A
        if A[1, 0] == 0 and A[2, 0] == 0:
            return True
        elif A[1, 0] % A[0, 0] == 0 and A[2, 0] % A[0, 0] == 0:
            self._first_finalize()
            self._Ps += self._L
            self._L = []
            return True
        else:
            return False

    def _first_one_loop(self):
        self._first_column()
        self._Ps += self._L
        self._L = []
        self._A = self._A.T
        self._first_column()
        self._Qs += self._L
        self._L = []
        self._A = self._A.T

    def _first_column(self):
        i = self._search_first_pivot()
        if i > 0:
            self._swap_rows(0, i)

        if self._A[1, 0] != 0:
            self._zero_first_column(1)
        if self._A[2, 0] != 0:
            self._zero_first_column(2)

    def _zero_first_column(self, j):
        if self._A[j, 0] < 0:
            self._flip_sign_row(j)
        A = self._A
        r, s, t = xgcd([A[0, 0], A[j, 0]])
        self._set_zero(0, j, A[0, 0], A[j, 0], r, s, t)

    def _search_first_pivot(self):
        A = self._A
        for i in range(3): # column index
            if A[i, 0] != 0:
                return i

    def _first_finalize(self):
        """Set zeros along the first colomn except for A[0, 0]

        This is possible only when A[1,0] and A[2,0] are dividable by A[0,0].

        """

        A = self._A
        L = np.eye(3, dtype='intc')
        L[1, 0] = -A[1, 0] // A[0, 0]
        L[2, 0] = -A[2, 0] // A[0, 0]
        self._L.append(L.copy())
        self._A = np.dot(L, self._A)

    def _second(self):
        """Find Smith normal form for Right-low 2x2 matrix"""

        self._second_one_loop()
        A = self._A
        if A[2, 1] == 0:
            return True
        elif A[2, 1] % A[1, 1] == 0:
            self._second_finalize()
            self._Ps += self._L
            self._L = []
            return True
        else:
            return False

    def _second_one_loop(self):
        self._second_column()
        self._Ps += self._L
        self._L = []
        self._A = self._A.T
        self._second_column()
        self._Qs += self._L
        self._L = []
        self._A = self._A.T

    def _second_column(self):
        """Right-low 2x2 matrix

        Assume elements in first row and column are all zero except for A[0,0].

        """

        if self._A[1, 1] == 0 and self._A[2, 1] != 0:
            self._swap_rows(1, 2)

        if self._A[2, 1] != 0:
            self._zero_second_column()

    def _zero_second_column(self):
        if self._A[2, 1] < 0:
            self._flip_sign_row(2)
        A = self._A
        r, s, t = xgcd([A[1, 1], A[2, 1]])
        self._set_zero(1, 2, A[1, 1], A[2, 1], r, s, t)

    def _second_finalize(self):
        """Set zero at A[2, 1]

        This is possible only when A[2,1] is dividable by A[1,1].

        """

        A = self._A
        L = np.eye(3, dtype='intc')
        L[2, 1] = -A[2, 1] // A[1, 1]
        self._L.append(L.copy())
        self._A = np.dot(L, self._A)

    def _swap_rows(self, i, j):
        """Swap i and j rows

        As the side effect, determinant flips.

        """

        L = np.eye(3, dtype='intc')
        L[i, i] = 0
        L[j, j] = 0
        L[i, j] = 1
        L[j, i] = 1
        self._L.append(L.copy())
        self._A = np.dot(L, self._A)

    def _flip_sign_row(self, i):
        """Multiply -1 for all elements in row"""

        L = np.eye(3, dtype='intc')
        L[i, i] = -1
        self._L.append(L.copy())
        self._A = np.dot(L, self._A)

    def _set_zero(self, i, j, a, b, r, s, t):
        """Let A[i, j] be zero based on Bezout's identity

           [ii ij]
           [ji jj] is a (k,k) minor of original 3x3 matrix.

        """

        L = np.eye(3, dtype='intc')
        L[i, i] = s
        L[i, j] = t
        L[j, i] = -b // r
        L[j, j] = a // r
        self._L.append(L.copy())
        self._A = np.dot(L, self._A)

def snf(mat):
    a = Snf(mat)
    a.run()
    return a.P, a.A, a.Q

def xgcd(vals):
    _xgcd = Xgcd(vals)
    return _xgcd.run()


class Xgcd(object):
    def __init__(self, vals):
        self._vals = numpy.array(vals, dtype='int')

    def run(self):
        r0, r1 = self._vals
        s0 = 1
        s1 = 0
        t0 = 0
        t1 = 1
        for i in range(1000):
            r0, r1, s0, s1, t0, t1 = self._step(r0, r1, s0, s1, t0, t1)
            if r1 == 0:
                break
        self._rst = numpy.array([r0, s0, t0], dtype='int')
        return self._rst

    def _step(self, r0, r1, s0, s1, t0, t1):
        q, m = divmod(r0, r1)
        r2 = m
        s2 = s0 - q * s1
        t2 = t0 - q * t1
        return r1, r2, s1, s2, t1, t2

    def __str__(self):
        v = self._vals
        r, s, t = self._rst
        return "%d = %d * (%d) + %d * (%d)" % (r, v[0], s, v[1], t)
