# -*- coding: utf-8 -*-
import numpy
import spglib
from pniggli import niggli_reduce

from itertools import product

from sagar.toolkit.mathtool import closest_pair
from sagar.element.base import periodic_table_dict, get_symbol, symbol2number

def car_to_frac(lattice, car_vec):
    inv_lattice = numpy.linalg.inv(lattice)
    return numpy.matmul(car_vec, inv_lattice)


def frac_to_car(lattice, frac_vec):
    return numpy.matmul(frac_vec, lattice)


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

    # For defining a immutable object
    def __hash__(self):
        pass

    def __eq__(self):
        pass

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
        positions = positions_1 - numpy.floor(positions_1)

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

        #m = numpy.amax(    mat)
        #_int_coor = numpy.array([i for i in product(range(m * 3), repeat=3)])
        # 适用于hnf

        # ---------------------选取一个大框包含目标框
        # import pdb; pdb.set_trace()
        conv = numpy.vstack((mat, numpy.array([0, 0, 0])))
        conv = numpy.vstack((conv, mat[0] + mat[1]))
        conv = numpy.vstack((conv, mat[0] + mat[2]))
        conv = numpy.vstack((conv, mat[1] + mat[2]))
        conv = numpy.vstack((conv, mat[0] + mat[1] + mat[2]))

        ma = [0, 0, 0]
        mi = [0, 0, 0]
        for i in range(0, 3):
            ma[i] = numpy.amax(conv[:, i])
            mi[i] = numpy.amin(conv[:, i])

        # import pdb; pdb.set_trace()
        _int_coor = numpy.array([j for j in product(range(
            mi[0] - 1, ma[0] + 2), range(mi[1] - 1, ma[1] + 2), range(mi[2] - 1, ma[2] + 2))])
        # from qiusb-------------

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

    def get_rotations_without_transitions(self, symprec=1e-5):
        """
        All rotations but remove inversion and its derive operation
        """
        rot_list = self.get_rotations(symprec)
        trans_list = self.get_pure_translations(symprec)
        rot_list_notrans = []
        for rot, trans in zip(rot_list, trans_list):
            if numpy.allclose(trans, [0, 0, 0]):
                rot_list_notrans.append(rot)
        return rot_list_notrans

    def get_pure_translations(self, symprec=1e-5):
        """
        get_pure_translations return pure translations.

        Info:
        For the supercell, translations and ratations can be seperated.
        total length of symmetry operations is len(rot)*len(trans)
        For the primitive cell, the translation is (0, 0, 0)
        """
        # trans_all = self.get_symmetry(symprec)['translations']
        # return numpy.unique(trans_all, axis=0)
        return self.get_symmetry(symprec)['translations']

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
        # 用下面这个原胞会改变坐标轴，no_idealize=True保持了坐标的方向。
        # lattice, positions, atoms = spglib.find_primitive(spg_cell, symprec)
        lattice, positions, atoms = spglib.standardize_cell(
            spg_cell, to_primitive=True, no_idealize=True, symprec=symprec)
        return self.__class__(lattice, positions, atoms)

    def _get_niggli_2D(self, vacc=16.0, eps=1e-5):
        # import pdb; pdb.set_trace()
        L = self.lattice[0:2, 0:2]
        reduced_L = niggli_reduce(L)
        # M = r_L . L^-1
        M = numpy.matmul(reduced_L, numpy.linalg.inv(L))
        M_3D = numpy.zeros((3, 3))
        M_3D[:2, :2] = M
        M_3D[2, 2] = 1
        M_3D = numpy.around(M_3D).astype(int)
        reduced_cell = self.extend(M_3D)

        lattice = reduced_cell.lattice
        lattice[2, 2] = vacc
        positions = reduced_cell.positions
        positions[:, 2] = 0.5
        atoms = reduced_cell.atoms

        return self.__class__(lattice, positions, atoms)

    def _get_niggli_3D(self, eps=1e-5):
        L = self.lattice
        reduced_L = niggli_reduce(L)
        # M = r_L . L^-1
        M = numpy.matmul(reduced_L, numpy.linalg.inv(L))
        reduced_cell = self.extend(M)

        return reduced_cell

    def get_refined_cell(self, symprec=1e-5):
        spg_cell = (self.lattice, self.positions, self.atoms)
        lattice, positions, atoms = spglib.refine_cell(spg_cell, symprec)
        return self.__class__(lattice, positions, atoms)

    def check(self, elements=None, limit=0.1, warn=False):
        """
        该方法用于自查对象中的位点是否过近
        若过近则抛出一个warning
        """
        if elements is None:
            origin_cell = self
        else:
            # 选取要check的元素，构建新的胞
            ele_num = [symbol2number(s) for s in elements]
            positions = []
            atoms = []
            for idx, e in enumerate(self.atoms):
                if e in ele_num:
                    positions.append(self.positions[idx])
                    atoms.append(self.atoms[idx])
            origin_cell = self.__class__(self.lattice, positions, atoms)

        mat = numpy.diag([3, 3, 3])
        periodic_cell = origin_cell.extend(mat)
        points = frac_to_car(periodic_cell.lattice, periodic_cell.positions)
        if closest_pair(points) < limit:
            if warn is True:
                import warnings
                warnings.warn("some atoms are too close(< {:f}), check cell".format(
                    limit), RuntimeWarning)
            return False
        else:
            return True

# 该类用于可以变换的结构
# TODO: 与上面的不可变的类同时继承于固定的基类
# TODO: 上面的类明确为不可变的类


class MutableCell(object):
    """
    MutableCell object represent a crystal structure.
    """

    def __init__(self, lattice, sites=None):
        """
        提供晶格和位点信息
        可以不提供位点信息，则初始化为空的晶胞
        在这个类中，site是一个list, list的第一个位置site[0]是坐标，第二个位置site[1]是原素类别
        """
        self._lattice = numpy.array(lattice).reshape((3, 3))
        if sites is None:
            self._sites = []
        else:
            self._sites = sites

    def to_cell(self):
        positions = []
        atoms = []
        for site in self._sites:
            positions.append(site[0])
            atoms.append(symbol2number(site[1]))
        cell = Cell(self._lattice,
                    numpy.array(positions).reshape((-1, 3)),
                    numpy.array(atoms))
        return cell

    @classmethod
    def from_cell(cls, imcell):
        lattice = imcell.lattice
        sites = []

        for p, e in zip(imcell.positions.tolist(), imcell.atoms.tolist()):
            sites.append([p, get_symbol(e)])

        return cls(lattice, sites)

    def add_site(self, site):
        """
        """
        self._sites.append(site)

    def set_site(self, idx, site):
        """
        """
        self._sites[idx] = site

    def remove_site(self, idx):
        """
        用set_site为空
        """
        pos = self._sites[idx][0]
        vacc_site = [pos, 'Vacc']
        self.set_site(idx, vacc_site)

    def rotate_site_by_z(self, idx, cc, degrees, radians=None):
        """
        需要确保是二维材料
        绕z旋转,则cc为(x,y),z方向确定不动.
        将idx的原子绕cc为圆心旋转angle度,angle单位是度,不是PI
        """
        radians = numpy.deg2rad(degrees)
        r_matrix = numpy.array([[numpy.cos(radians), -numpy.sin(radians)],
                                [numpy.sin(radians), numpy.cos(radians)]])
        x, y, z = frac_to_car(self._lattice, self._sites[idx][0])
        e = self._sites[idx][1]
        x0, y0, _ = frac_to_car(self._lattice, cc)
        # 原来的向量
        vx, vy = x-x0, y-y0
        # 旋转后的向量
        new_vx, new_vy = numpy.matmul([vx, vy], r_matrix)
        new_x, new_y = new_vx+x0, new_vy+y0
        new_site = [car_to_frac(self._lattice, (new_x, new_y, z)), e]
        self.set_site(idx, new_site)

    def get_site(self, idx):
        return self._sites[idx]

    def get_sites(self, start=None, stop=None, step=1):
        """
        使用slice(start, stop, step)获取多个原子.
        """
        idxs = slice(start, stop, step)
        return self._sites[idxs]

    def get_car_site(self, idx):
        fs = self.get_site(idx)
        car = frac_to_car(self._lattice, fs[0])

        return [tuple(car), fs[1]]

    def __repr__(self):
        def _repr(number):
            return "{:9.6f}".format(number)

        lattice = self._lattice
        out_latt = ["Lattice:",
                    "   a: " + ' '.join(map(_repr, lattice[0])),
                    "   b: " + ' '.join(map(_repr, lattice[1])),
                    "   c: " + ' '.join(map(_repr, lattice[2]))]
        out_pos = []
        out_pos.append("Sites:")
        for idx, s in enumerate(self._sites):
            idx = '[' + str(idx) + ']:'
            o = idx + ' '.join(map(_repr, s[0])) + ' ' + s[1]
            out_pos.append(o)

        outs = out_latt + out_pos
        return "\n".join(outs)
