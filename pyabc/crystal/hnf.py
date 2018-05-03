import numpy

from math import sqrt
from itertools import product

from progressbar import ProgressBar


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
    # TODO: first make sure pcell is primitive
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
    """
    for rot in rot_list:
        m = numpy.matmul(
            numpy.matmul(hnf_x, numpy.linalg.inv(rot.T)),
            numpy.linalg.inv(hnf_y))
        if numpy.allclose(numpy.mod(m, 1), numpy.zeros_like(m), atol=prec):
            return True

    return False
