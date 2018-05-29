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

def symbol2number(symbol):
    return periodic_table_dict[symbol]

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
        # TODO: magmoms init setting
        self._lattice = numpy.array(lattice).reshape((3, 3))

        self._atom_numbers = len(atoms)
        # import pdb; pdb.set_trace()
        self._positions = numpy.array(positions).reshape((-1, 3))
        moded = numpy.ones_like(self._positions, dtype='intc')
        self._positions = numpy.mod(self._positions, moded)
        if self._positions.shape[0] != self._atom_numbers:
            raise ValueError("When init Cell, number of atoms not equal to "
                             "number of positions.\n"
                             "CHECK YOUR INPUT!")

        if isinstance(atoms, numpy.ndarray):
            atoms = atoms.tolist()
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
            else:

                a.append(round(s))

        self._atoms = numpy.array(a, dtype='intc')
        # TODO: initial with Cartesian coor

    @property
    def lattice(self):
        return self._lattice

    @property
    def positions(self):
        return self._positions

    @property
    def atoms(self):
        return self._atoms

    def __repr__(self):
        def _repr(number):
            return "{:9.6f}".format(number)
        out_latt = ["Lattice:",
                    "   a: " + ' '.join(map(_repr, self._lattice[0])),
                    "   b: " + ' '.join(map(_repr, self._lattice[1])),
                    "   c: " + ' '.join(map(_repr, self._lattice[2]))]
        sites = zip(self._positions, self._atoms)
        out_pos = []
        out_pos.append("Sites:")
        for s in sites:
            o = ' '.join(map(_repr, s[0])) + ' ' + get_symbol(s[1])
            out_pos.append(o)

        outs = out_latt + out_pos
        return "\n".join(outs)

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
        #       2. 针对非对角矩阵，是一样适用？ (DONE)
        # TODO: 不必是hnf矩阵 (DONE)
        # TODO: mat必须是一个整数矩阵，做一个判断，给出异常
        lattice = numpy.matmul(mat, self._lattice)

        smallest_cell = numpy.matmul(self._positions, numpy.linalg.inv(mat))
        grids = self._get_mat_frac(mat)
        list_positions = [i for i in map(
            lambda x: x + grids, list(smallest_cell))]
        positions_1 = numpy.concatenate(list_positions, axis=0)
        positions = positions_1-numpy.floor(positions_1)

        n = abs(numpy.linalg.det(mat))
        atoms = numpy.repeat(self._atoms, int(round(n)))

        return self.__class__(lattice, positions, atoms)

    def _get_mat_frac(self, mat):
        """
        When giving mat -- a 3x3 matrix,
        export a numpy.array represent the
        grid points between 0~1.

        Used in producing the new positions extended by a matrix
        """
        prec = 1e-5
        #m = numpy.amax(mat)
        #_int_coor = numpy.array([i for i in product(range(m * 3), repeat=3)])
        ##适用于hnf


        #---------------------选取一个大框包含目标框
        conv =  numpy.row_stack((mat, numpy.matrix([0, 0, 0])))
        conv =  numpy.row_stack((conv, numpy.matrix(mat[0]+mat[1])))
        conv =  numpy.row_stack((conv, mat[0]+mat[2]))
        conv =  numpy.row_stack((conv, mat[1]+mat[2]))
        conv =  numpy.row_stack((conv, mat[0]+mat[1]+mat[2]))

        ma = [0, 0, 0]
        mi = [0, 0, 0]
        for i in range(0,3):
            ma[i] = numpy.amax(conv[:,i])
            mi[i] = numpy.amin(conv[:,i])
        _int_coor =  numpy.array([j for j in product(range(mi[0]-1, ma[0]+2), range(mi[1]-1, ma[1]+2), range(mi[2]-1, ma[2]+2))])
        ##from qiusb-------------

        _all_frac = numpy.matmul(_int_coor, numpy.linalg.inv(mat))

        is_incell = numpy.all(
            ((_all_frac > -prec) & (_all_frac < 1 - prec)), axis=1)

        return _all_frac[numpy.where(is_incell)[0]]

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

    def get_rotations_without_inversion(self, symprec=1e-5):
        """
        All rotations but remove inversion and its derive operation
        """
        rot_list = self.get_rotations(symprec)
        rot_list_noinv = []
        for rot in rot_list:
            if self._not_contain(rot_list_noinv, rot):
                rot_list_noinv.append(rot)
        return rot_list_noinv

    def get_pure_translations(self, symprec=1e-5):
        """
        get_pure_translations return pure translations.

        Info:
        For the supercell, translations and ratations can be seperated.
            total length of symmetry operations is len(rot)*len(trans)
        For the primitive cell, the translation is (0, 0, 0)
        """
        trans_all = self.get_symmetry(symprec)['translations']
        return numpy.unique(trans_all, axis=0)

    def _not_contain(self, rot_list_noinv, rot):
        sym_inv = numpy.array([-1, 0, 0,
                               0, -1, 0,
                               0, 0, -1]).reshape((3, 3))
        for r in rot_list_noinv:
            rot_inv = numpy.matmul(rot, sym_inv)
            if numpy.allclose(r, rot) or numpy.allclose(r, rot_inv):
                return False
        return True

    def get_spacegroup(self, symprec=1e-5):
        """
        dependent on spglib https://atztogo.github.io/spglib/
        """
        return spglib.get_spacegroup((self._lattice, self._positions, self._atoms), symprec)

    def is_primitive(self, symprec=1e-5):
        """
        is_primitive_cell decide if a cell is primitive

        parameters:

        cell: Cell object
        prec: float, the precision to judge, default=1e-5

        return: bool
        """
        natoms = len(self.atoms)
        spg_cell = (self.lattice, self.positions, self.atoms)
        pnatoms = len(spglib.find_primitive(spg_cell, symprec)[2])
        return natoms == pnatoms

    def get_primitive_cell(self, symprec=1e-5):
        spg_cell = (self.lattice, self.positions, self.atoms)
        lattice, positions, atoms = spglib.find_primitive(spg_cell, symprec)
        return self.__class__(lattice, positions, atoms)

    def get_refine_cell(self, symprec=1e-5):
        spg_cell = (self.lattice, self.positions, self.atoms)
        lattice, positions, atoms = spglib.refine_cell(spg_cell, symprec)
        return self.__class__(lattice, positions, atoms)
