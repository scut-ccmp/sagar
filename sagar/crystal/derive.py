# -*- coding: utf-8 -*-
import numpy

from sagar.crystal.utils import non_dup_hnfs, snf
from sagar.toolkit.mathtool import is_int_np_array, refine_positions
from sagar.toolkit.derivetool import remove_redundant

from sagar.crystal.structure import Cell


def cells_nonredundant(pcell, volume=1, dimension=3, symprec=1e-5, comprec=1e-5):
    """
    cells_nonredundant return all non duplicated hnf extend cells.

    parameters:

    pcell: Cell object, The primitive cell to be extended
    volume: int, Extend to how large supercellself, default=1
    symprec: int, symmetry precision
    When finding duplicated hnfs the precesion, default=1e-5
    comprec: float, compare precision
    When finding the rotations symmetry of primitive cell, defalut=1e-5

    yield:
    A list of Cell objects.
    """
    # return [pcell.extend(hnf) for hnf in non_dup_hnfs(pcell, volume, dimension, symprec, comprec)]
    for hnf in non_dup_hnfs(pcell, volume, dimension, symprec, comprec):
        cell = pcell.extend(hnf)
        if dimension == 2:
            cell = cell._get_niggli_2D(vacc=16, eps=comprec)
        if dimension == 3:
            cell = cell.get_refined_cell(symprec=symprec)
        yield cell


class PermutationGroup(object):
    """
    所有的对称操作都是以置换矩阵的形式，作用在一个元素排列上。
    """

    def __init__(self, pcell, mat):
        if not isinstance(pcell, Cell):
            raise TypeError(
                "want sagar.crystal.structure.Cell, got {:}".format(type(pcell)))
        self._pcell = pcell
        self._mat = mat
        self._snf, L, R = snf(mat)
        self._volume = numpy.diagonal(self._snf).prod()
        self._nsites = len(pcell.atoms)  # 最小原胞中原子个数 如：hcp为2

    def get_symmetry_perms(self, symprec=1e-5):
        # Q??: use whose rotations?????
        supercell = self._pcell.extend(self._mat)
        # import pdb; pdb.set_trace()
        # 用超胞的旋转对称才是合理的
        syms = supercell.get_symmetry(symprec)
        arr_rots, arr_trans = syms['rotations'], syms['translations']

        result = numpy.zeros(
            (len(arr_rots), self._nsites * self._volume), dtype='intc')

        origin_positions = supercell.positions
        origin_positions = refine_positions(origin_positions)
        for i, (rot, trans) in enumerate(zip(arr_rots, arr_trans)):
            new_positions = numpy.matmul(origin_positions, rot.T) + trans
            moded = numpy.ones_like(new_positions, dtype='intc')
            new_positions = numpy.mod(new_positions, moded)
            new_positions = refine_positions(new_positions, atol=symprec)
            # 寻找置换矩阵
            for j, row in enumerate(origin_positions):
                row = refine_positions(row, atol=symprec)
                idx = numpy.where(
                    (numpy.isclose(row, new_positions, atol=symprec)).all(axis=1))[0]
                result[i, j] = idx

        result = numpy.unique(result, axis=0)
        return result

# 在较大体积的构型中，可能有结构已经在小的体积中出现过。
# 那么是否要移除这些已经出现过的结构呢？
#
# 暂时给出的答案是：若要目的是产生基态相图，则出现过的无需重复计算。
#     若目的是计算结构能量做统计平均，则在每个体积倍数下需要得到简并度，以及总体数目，
#     所以该体积倍数下的所有构型都是需要的，包括完全不取代的构型。
#     TODO: 在以后的发展中，可以用到结构对比的技术，以允许统计的情况下无需重复计算。
#
# 该函数中没有comprec参数，是因为构型(configurations)之间的对比
# 通过原子排列的对比实现，没有坐标的比较。


class ConfigurationGenerator(object):
    """
    pcell: Cell object, The primitive cell to be extended
    """

    def __init__(self, cell, symprec=1e-5):
        self._cell = cell
        if not cell.is_primitive(symprec):
            self._pcell = cell.get_primitive_cell(symprec)
        else:
            self._pcell = cell

    def cons_max_volume(self, sites, max_volume, min_volume=1, dimension=3, symprec=1e-5):
        """
        parameters:

        pcell: Cell object, The primitive cell to be extended
        sites: disorderd sites infomation.

        yield:

        Cell object, a list of non-redundant configurations to max volume supercells.

        """
        # 该函数产生所有构型用于确定基态相图
        for volume in range(min_volume, max_volume + 1):
            hnfs = non_dup_hnfs(self._pcell, volume, dimension, symprec)
            dict_trans = {}  # 记录已经产生过的snf，相同snf的平移操作相同。
            for h in hnfs:
                hfpg = PermutationGroup(self._pcell, h)
                # 此处的平移操作不必每次重新计算，因为相同snf平移操作相同
                # 可以用字典来标记查询。若没有这样的操作，那么就没有snf带来的效率提升。
                # quotient = hfpg.get_quotient()
                # if not quotient in dict_trans:
                #     dict_trans[quotient] = hfpg.get_pure_translations(symprec)
                # trans = dict_trans[quotient]
                # rots = hfpg.get_pure_rotations(symprec)
                perms = hfpg.get_symmetry_perms(symprec)

                supercell = self._pcell.extend(h)
                _sites = numpy.repeat(sites, volume, axis=0)

                for mol, _ in remove_redundant(supercell.positions, _sites, perms):
                    c = Cell(supercell.lattice, mol[0], mol[1])
                    if c.is_primitive(symprec):
                        yield c

    # 特定体积胞
    def cons_specific_volume(self, sites, volume=2, e_num=None, dimension=3, symprec=1e-5):
        """
        parameters:

        sites: 2D list, disorderd sites infomation.
        volume: int, certain volume with subperiodic.
        symprec: float, precision for symmetry find.

        yield:

        a tuple
        tuple[0]: Cell object, a list of non-redundant configurations of certain volume supercell.
        tuple[1]: int object, degeneracy of the configuration in all configurations of this volume.
        """
        # TODO: if it is a supercell as input: get Error
        # 该函数产生特定体积下所有构型（包括超胞）和简并度，用于统计平均
        hnfs = non_dup_hnfs(self._pcell, volume, dimension, symprec)
        dict_trans = {}
        for h in hnfs:
            hfpg = PermutationGroup(self._pcell, h)
            # 此处的平移操作不必每次重新计算，因为相同snf平移操作相同
            # 可以用字典来标记查询。若没有这样的操作，那么就没有snf带来的效率提升。
            # For hnf with same snf, translations are same.
            # quotient = hfpg.get_quotient()
            # if not quotient in dict_trans:
            #     dict_trans[quotient] = hfpg.get_pure_translations(symprec)
            # trans = dict_trans[quotient]
            # rots = hfpg.get_pure_rotations(symprec)
            perms = hfpg.get_symmetry_perms(symprec)

            supercell = self._pcell.extend(h)
            _sites = numpy.repeat(sites, volume, axis=0)

            for mol, d in remove_redundant(supercell.positions, _sites, perms, e_num):
                c = Cell(supercell.lattice, mol[0], mol[1])
                yield (c, d)

    def cons_specific_cell(self, sites, e_num=None, symprec=1e-5):
        """
        cons_specific_cell_and_c generate configurations of specific cell
        and specific concentration.

        parameters:
        sites: list of (lists or tuples), represent element disorder of each sites
        e_num: tuple, number of atoms in disorderd sites.

        e_num 特指无序位点的浓度，也就是原子比，而不是整体构型的元素原子比。有可能其他位点存在相同构型。
        !!限制，无序位点的组成必须是相同的，而上面几个函数的无序位点是可以不同的。!!
        """
        lat_cell = self._cell.lattice
        lat_pcell = self._pcell.lattice
        mat = numpy.matmul(lat_cell, numpy.linalg.inv(lat_pcell))
        if is_int_np_array(mat):
            mat = numpy.around(mat).astype('intc')
        else:
            print("cell:\n", lat_cell)
            print("primitive cell:\n", lat_pcell)
            raise ValueError(
                "cell lattice and its primitive cell lattice not convertable")
        hfpg = PermutationGroup(self._pcell, mat)

        perms = hfpg.get_symmetry_perms(symprec)

        supercell = self._pcell.extend(mat)
        for mol, d in remove_redundant(supercell.positions, sites, perms, e_num=e_num):
            c = Cell(supercell.lattice, mol[0], mol[1])
            yield (c, d)
