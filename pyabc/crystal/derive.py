import numpy
from itertools import product

from pyabc.crystal.utils import non_dup_hnfs, is_int_np_array, binomialCoeff
from pyabc.crystal.utils import snf, refine_positions

from pyabc.crystal.structure import Cell


def cells_nonredundant(pcell, volume=1, symprec=1e-5, comprec=1e-5):
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
    return [pcell.extend(hnf) for hnf in non_dup_hnfs(pcell, volume, symprec, comprec)]
    # for hnf in non_dup_hnfs(pcell, volume, symprec, comprec):
    #     yield pcell.extend(hnf)

class PermutationGroup(object):
    """
    所有的对称操作都是以置换矩阵的形式，作用在一个元素排列上。
    """

    def __init__(self, pcell, mat):
        if not isinstance(pcell, Cell):
            raise TypeError(
                "want pyabc.crystal.structure.Cell, got {:}".format(type(cell)))
        self._pcell = pcell
        self._mat = mat
        self._snf, L, R = snf(mat)
        self._volume = numpy.diagonal(self._snf).prod()
        self._nsites = len(pcell.atoms)  # 最小原胞中原子个数 如：hcp为2

    def get_symmetry_perms(self, symprec=1e-5):
        # Q??: use whose rotations?????
        supercell = self._pcell.extend(self._mat)
        # 用超胞的旋转对称才是合理的
        # arr_rots = supercell.get_rotations(symprec)[:]  # 第一个是单位矩阵
        arr_rots = supercell.get_rotations(symprec)[:]  # 第一个是单位矩阵
        arr_trans = supercell.get_pure_translations(symprec)[:]  # 第一个是单位矩阵
        result = numpy.zeros(
            (len(arr_rots), self._nsites * self._volume), dtype='intc')

        origin_positions = supercell.positions
        origin_positions = refine_positions(origin_positions)
        for i, (rot, trans) in enumerate(zip(arr_rots, arr_trans)):
            new_positions = numpy.matmul(origin_positions, rot.T) + trans
            moded = numpy.ones_like(new_positions, dtype='intc')
            new_positions = numpy.mod(new_positions, moded)
            new_positions = refine_positions(new_positions)
            # 寻找置换矩阵
            for j, row in enumerate(origin_positions):
                row = refine_positions(row)
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

    def cons_max_volume(self, sites, max_volume, min_volume=1, symprec=1e-5):
        """
        parameters:

        pcell: Cell object, The primitive cell to be extended
        sites: disorderd sites infomation.

        yield:

        Cell object, a list of non-redundant configurations to max volume supercells.

        """
        # 该函数产生所有构型用于确定基态相图
        for volume in range(min_volume, max_volume + 1):
            hnfs = non_dup_hnfs(self._pcell, volume, symprec)
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

                for c, _ in self._remove_redundant(supercell, sites, perms, volume, remove_super=True):
                    yield c

    # 特定体积胞
    def cons_specific_volume(self, sites, volume=2, symprec=1e-5):
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
        # 该函数产生特定体积下所有构型（包括超胞）和简并度，用于统计平均
        hnfs = non_dup_hnfs(self._pcell, volume, symprec)
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

            for c, d in self._remove_redundant(supercell, sites, perms, volume):
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

        for c, d in self._remove_redundant(self._cell, sites, perms, e_num=e_num):
            yield (c, d)

    def _remove_redundant(self, cell, sites, perms, volume=1, e_num=None, remove_super=False):
        # perms = self._get_perms_from_rots_and_trans(rots, trans)
        # TODO: 加入一个机制，来清晰的设定位点上无序的状态
        arg_sites = [len(i) for i in sites]
        arg_sites = numpy.repeat(arg_sites, volume)
        # redundant configurations do not want see again
        # 用于记录在操作作用后已经存在的构型排列，而无序每次都再次对每个结构作用所有操作
        redundant = set()

        # deg_total = 0
        # loop over configurations
        for atoms_mark in _atoms_gen(arg_sites, e_num):
            arr_atoms_mark = numpy.array(atoms_mark)
            ahash = _hash_atoms(atoms_mark)

            if remove_super:
                flag = (ahash in redundant) or self._is_super(
                    cell, arr_atoms_mark)
            else:
                flag = ahash in redundant

            if flag:
                continue
            else:
                list_all_transmuted = []
                for p in perms:
                    atoms_transmuted = arr_atoms_mark[p]
                    redundant.add(_hash_atoms(atoms_transmuted))
                    # degeneracy
                    list_all_transmuted.append(atoms_transmuted)

                arr_all_transmuted = numpy.array(list_all_transmuted)
                deg = numpy.unique(arr_all_transmuted, axis=0).shape[0]

                atoms = self._mark_to_atoms(arr_atoms_mark, sites)
                # print(str(atoms) + '  ' + str(deg))

                c = Cell(cell.lattice, cell.positions, atoms)
                yield (c, deg)

        #     deg_total += deg
        # print(deg_total)

    def _get_perms_from_rots_and_trans(self, rots, trans):
        nrot = rots.shape[0]
        ntran = trans.shape[0]
        size = nrot * ntran
        perms = numpy.zeros((size, trans.shape[1]), dtype='int')
        idx = 0
        for r in rots:
            for t in trans:
                p = r[t]
                perms[idx] = p
                idx += 1
        return perms

    def _is_super(self, cell, arr_atoms_mark):
        newcell = Cell(cell.lattice, cell.positions, arr_atoms_mark)
        return not newcell.is_primitive()

    def _mark_to_atoms(self, arr_mark, sites):
        num_of_site_groups = len(sites)
        arr_atoms = arr_mark.reshape(num_of_site_groups, -1)
        # import pdb; pdb.set_trace()
        atoms = numpy.zeros_like(arr_atoms)
        for i, row in enumerate(arr_atoms):
            for j, v in enumerate(row):
                atoms[i][j] = sites[i][v]

        return atoms.flatten().tolist()


def _atoms_gen(args, e_num=None):
    """
    parameter:
    args: list, represent number of atoms of each site.
    e_num: None or tuple, concentration of element in disorderd sites,
            default=None, for all configurations

    给定每个位点原子无序的数目，产生所有可能的原子排列，共k^n种。
    TODO: 在这里加入浓度比？
    """
    if e_num is None:
        p = []
        for i in args:
            p.append(range(i))
        return product(*p)
    else:
        arr_arrange = _serial_int_to_arrangement(e_num)
        # import pdb; pdb.set_trace()
        for col, n in enumerate(args):
            if n == 1:
                arr_arrange = numpy.insert(arr_arrange, col, 0, axis=1)
        return arr_arrange


def _hash_atoms(atoms):
    """
    给出一个atoms排列，返回一个string代表该排列。可用于hash表。
    """
    return ''.join(str(i) for i in atoms)


def _serial_int_to_arrangement(e_num):
    """
    Algorithm From:
    Hart, G. L. W., Nelson, L. J., & Forcade, R. W. (2012).
    Generating derivative structures at a fixed concentration, 59, 101–107.

    Corrections:
    Some error in paper: fig. 5 --- loop over site: should be t = t - 1
    """
    slots_total = _slots_total = sum(e_num)
    comb = []
    max = 1
    for v in e_num:
        bino = binomialCoeff(_slots_total, v)
        comb.append(bino)
        max *= bino
        _slots_total -= v

    max = int(max)
    arr_arrangement = numpy.full((max, slots_total), -1)

    for i in range(max):
        open_slots = slots_total
        for e in range(len(e_num)):
            y, x = divmod(i, comb[e])
            a = e_num[e]
            m = open_slots

            for idx, value in enumerate(arr_arrangement[i]):
                if value != -1:
                    # 该slot已经放置原子
                    continue
                else:
                    b = binomialCoeff(m - 1, a - 1)
                    if b <= x:
                        x -= b
                    else:
                        a -= 1
                        arr_arrangement[i][idx] = e
                    m -= 1

    return arr_arrangement
