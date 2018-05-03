import numpy

from itertools import product
import spglib


class Cell(object):
    """
    Cell object represent a crystal structure.

    parameters:

    lattice: 3x3 1D-list, lattice of cell

    positions: n tuples(x,y,z) in fraction. 可多态初始化

    atoms: list of numbers, represent atom in periodic table.
    """

    def __init__(self, lattice, positions, atoms):
        self._lattice = numpy.array(lattice).reshape((3, 3))
        self._positions = numpy.array(positions)
        self._atoms = numpy.array(atoms)
        self._numbers = len(atoms)

    def extend(self, mat):
        lattice = numpy.matmul(mat, self._lattice)

        smallest_cell = numpy.matmul(self._positions, numpy.linalg.inv(mat))
        grids = _get_mat_frac(mat)
        list_positions = [i for i in map(
            lambda x: x + grids, list(smallest_cell))]
        positions = numpy.concatenate(list_positions, axis=0)

        n = mat.diagonal().prod()
        atoms = numpy.repeat(self._numbers, n)

        return self.__class__(lattice, positions, atoms)

    def get_symmetry(self, symprec=1e-5):
        """
        dependent on spglib https://atztogo.github.io/spglib/
        """
        return spglib.get_symmetry((self._lattice, self._positions, self._atoms), symprec)

    def get_rotations(self, symprec=1e-5):
        """
        dependent on spglib https://atztogo.github.io/spglib/
        """
        return spglib.get_symmetry((self._lattice, self._positions, self._atoms), symprec)['rotations']

    def get_spacegroup(self, symprec=1e-5):
        """
        dependent on spglib https://atztogo.github.io/spglib/
        """
        return spglib.get_spacegroup((self._lattice, self._positions, self._atoms), symprec)


def _get_mat_frac(mat):
    """
    When giving mat -- a 3x3 matrix,
    export a numpy.array represent the
    grid points between 0~1.

    Used in producing the new positions extended by a matrix
    """
    prec = 1e-5
    m = numpy.amax(mat)
    _int_coor = numpy.array([i for i in product(range(m * 3), repeat=3)])
    _all_frac = numpy.matmul(_int_coor, numpy.linalg.inv(mat))

    is_incell = numpy.all(
        ((_all_frac > -prec) & (_all_frac < 1 - prec)), axis=1)

    return _all_frac[numpy.where(is_incell)[0]]
