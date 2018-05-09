import numpy

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
        is_array_int = numpy.all(numpy.isclose(m, numpy.around(m), atol=1e-3))
        if is_array_int:
            return True

    return False


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
    opL = numpy.eye(3, dtype='int')
    opR = numpy.eye(3, dtype='int')
    while not _is_incremental_diag(mat):
        pivot = _search_first_pivot(mat)
        if pivot > 0:
            mat, op = _swap_rows(mat, 0, pivot)
            opL = numpy.matmul(op, opL)

        if mat[1, 0] % mat[0, 0] == 0 and mat[2, 0] % mat[0, 0] == 0:
            mat, op = _first_exact_division(mat)
            opL = numpy.matmul(op, opL)

        mat, op = _zero_first_column(mat)
        opL = numpy.matmul(op, opL)

        mat, op = _zero_first_row(mat)
        opR = numpy.matmul(opR, op)

        # lower 2x2 matrix
        if mat[1, 1] == 0 and mat[2, 1] != 0:
            mat, op = _swap_rows(mat, 1, 2)
            opL = numpy.matmul(op, opL)

        if mat[2, 1] % mat[1, 1] == 0:
            mat, op = _second_exact_division(mat)
            opL = numpy.matmul(op, opL)

        mat, op = _zero_second_column(mat)
        opL = numpy.matmul(op, opL)

        mat, op = _zero_second_row(mat)
        opR = numpy.matmul(opR, op)

        if mat[2, 2] < 0:
            op = numpy.array([1, 0, 0,
                              0, 1, 0,
                              0, 0, -1]).reshape((3, 3))
            opL = numpy.matmul(op, opL)
            mat = numpy.matmul(op, mat)

        # sort the diag and get the opL matrix
        mat_flat = numpy.diagonal(mat)
        idx = numpy.argsort(mat_flat)
        op = numpy.eye(3, dtype='int')[idx]
        mat_tmp = numpy.matmul(op, mat)
        opL = numpy.matmul(op, opL)

        mat = numpy.matmul(mat_tmp, op.T)
        opR = numpy.matmul(opR, op.T)

    return opL, mat, opR


def _is_diag(mat):
    return numpy.all(mat == numpy.diag(numpy.diagonal(mat)))


def _is_incremental_diag(mat):
    is_diag = numpy.all(mat == numpy.diag(numpy.diagonal(mat)))
    list_mat_flat = numpy.diagonal(mat).tolist()
    is_incremental = sorted(list_mat_flat) == list_mat_flat
    return is_diag and is_incremental


def _search_first_pivot(mat):
    for i in range(3):
        if mat[i, 0] != 0:
            return i


def _swap_rows(mat, i, j):
    """
    return:
    mat: new row permuted mat
    op: is a row permute matrix
    """
    op = numpy.eye(3, dtype='int')
    op[i, i], op[j, j] = 0, 0
    op[i, j], op[j, i] = 1, 1
    mat = numpy.matmul(op, mat)
    return mat, op


def _zero_first_column(mat):
    op_return = numpy.eye(3, dtype='int')
    for i in [1, 2]:
        if mat[i, 0] != 0:
            mat, op = _zero_first_ele_in_row_i(mat, i)
            op_return = numpy.matmul(op, op_return)
    return mat, op_return


def _zero_first_ele_in_row_i(mat, i):
    op_return = numpy.eye(3, dtype='int')
    if mat[i, 0] < 0:
        mat, op = _flip_sign_row(mat, i)
        op_return = numpy.matmul(op, op_return)
    r, s, t = extended_gcd(mat[0, 0], mat[i, 0])
    mat, op = _set_zero(mat, 0, i, mat[0, 0], mat[i, 0], r, s, t)
    op_return = numpy.matmul(op, op_return)
    return mat, op_return


def _first_exact_division(mat):
    op = numpy.eye(3, dtype='int')
    op[1, 0] = -mat[1, 0] // mat[0, 0]
    op[2, 0] = -mat[2, 0] // mat[0, 0]
    mat = numpy.matmul(op, mat)
    return mat, op


def _zero_first_row(mat):
    mat, op = _zero_first_column(mat.T)
    return mat.T, op.T


def _zero_second_column(mat):
    op_return = numpy.eye(3, dtype='int')
    if mat[2, 1] < 0:
        mat, op = _flip_sign_row(mat, 2)
        op_return = numpy.matmul(op, op_return)
    r, s, t = extended_gcd(mat[1, 1], mat[2, 1])
    mat, op = _set_zero(mat, 1, 2, mat[1, 1], mat[2, 1], r, s, t)
    op_return = numpy.matmul(op, op_return)
    return mat, op_return


def _second_exact_division(mat):
    op = numpy.eye(3, dtype='int')
    op[2, 1] = -mat[2, 1] // mat[1, 1]
    mat = numpy.matmul(op, mat)
    return mat, op


def _zero_second_row(mat):
    mat, op = _zero_second_column(mat.T)
    return mat.T, op.T


def _flip_sign_row(mat, i):
    """
    mul -1 for i row
    """
    op = numpy.eye(3, dtype='int')
    op[i, i] = -1
    mat = numpy.matmul(op, mat)
    return mat, op


def _set_zero(mat, i, j, aa, bb, r, s, t):
    """
    Based on Bezout's identity

    Let mat[i, j] be zero
    """
    op = numpy.eye(3, dtype='int')
    op[i, i] = s
    op[i, j] = t
    op[j, i] = -bb // r
    op[j, j] = aa // r
    mat = numpy.matmul(op, mat)
    return mat, op
