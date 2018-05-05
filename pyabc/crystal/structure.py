import numpy

from itertools import product
import spglib

"""
periodic_table_dict
"""
periodic_table_dict = {'Vacc': 0,
                       'H': 1, 'He': 2,
                       'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
                       'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, }


def get_symbol(atom):
    """
    get_symbol return symbol of atomic number

    parameter:

    atom: int, atomic number.
    """
    for key, value in periodic_table_dict.items():
        if atom == value:
            return str(key)
    return "NaN_x"


class Cell(object):
    """
    Cell object represent a crystal structure.

    parameters:

    lattice: 3x3 1D-list, lattice of cell
    positions: n tuples(x,y,z) in fraction. 可多态初始化
    atoms: list of atoms, can be atomic number (int), can be atomic symbol (string),
           represent atom in periodic table. (用全称初始化？？谁TM会这么用？)
    """

    def __init__(self, lattice, positions, atoms):
        self._lattice = numpy.array(lattice).reshape((3, 3))

        self._atom_numbers = len(atoms)
        self._positions = numpy.array(positions).reshape((-1, 3))
        if self._positions.shape[0] != self._atom_numbers:
            raise ValueError("When init Cell, number of atoms not equal to "
                             "number of positions.\n"
                             "CHECK YOUR INPUT!")

        a = []
        for s in atoms:
            if isinstance(s, str):
                if s in periodic_table_dict:
                    a.append(periodic_table_dict[s])
                elif len(s.split("_")) == 2 and s.split("_")[0] == "NaN":
                    i = int(s.split("_")[1])
                    a.append(1000 + i)
                else:
                    raise ValueError("Unkown atom symbols {:}".format(s))
            elif isinstance(s, int):
                a.append(s)

        self._atoms = numpy.array(a)
        # TODO: initial with Cartesian coor

        self.lattice = self._lattice
        self.positions = self._positions
        self.atoms = self._atoms

    def extend(self, mat):
        """
        extend transform a cell by mat

        parameters:
        mat: a 3x3 numpy.array

        return:
        A new Cell object.
        """
        # TODO: now extend is proved right only for hnf matrixself.
        #       1. 我们是否需要把旋转合并进来？
        #       2. 针对非对角矩阵，是一样适用？
        lattice = numpy.matmul(mat, self._lattice)

        smallest_cell = numpy.matmul(self._positions, numpy.linalg.inv(mat))
        grids = _get_mat_frac(mat)
        list_positions = [i for i in map(
            lambda x: x + grids, list(smallest_cell))]
        positions = numpy.concatenate(list_positions, axis=0)

        n = mat.diagonal().prod()
        atoms = numpy.repeat(self._atoms, n)

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


def is_primitive_cell(cell, prec=1e-5):
    """
    is_primitive_cell decide if a cell is primitive

    parameters:

    cell: Cell object
    prec: float, the precision to judge, default=1e-5

    return: bool
    """
    natoms = len(cell.atoms)
    spg_cell = (cell.lattice, cell.positions, cell.atoms)
    pnatoms = len(spglib.find_primitive(spg_cell, prec)[2])
    return natoms == pnatoms
