import numpy

from pyabc.crystal.utils import non_dup_hnfs, atoms_gen, hash_atoms
from pyabc.crystal.utils import HartForcadePermutationGroup as HFPG

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
    # return [pcell.extend(hnf) for hnf in non_dup_hnfs(pcell, volume, symprec, comprec)]
    for hnf in non_dup_hnfs(pcell, volume, symprec, comprec):
        yield pcell.extend(hnf)


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
        for v in range(min_volume, max_volume + 1):
            hnfs = non_dup_hnfs(self._pcell, v, symprec)
            dict_trans = {}  # 记录已经产生过的snf，相同snf的平移操作相同。
            for h in hnfs:
                hfpg = HFPG(self._pcell, h)
                # 此处的平移操作不必每次重新计算，因为相同snf平移操作相同
                # 可以用字典来标记查询。若没有这样的操作，那么就没有snf带来的效率提升。
                quotient = hfpg.get_quotient()
                if not quotient in dict_trans:
                    dict_trans[quotient] = hfpg.get_pure_translations(symprec)
                trans = dict_trans[quotient]

                supercell = hfpg.get_supercell()

                # 产生所有可能操作的置换操作
                # perms = hfpg.get_symmetry()
                rots = hfpg.get_pure_rotations(symprec)
                perms = self._get_perms_from_rots_and_trans(rots, trans)

                # TODO: 加入一个机制，来清晰的设定位点上无序的状态
                arg_sites_pcell = [len(i) for i in sites]
                arg_sites = numpy.repeat(arg_sites_pcell, v)

                # redundant configurations do not want see again
                # 用于记录在操作作用后已经存在的构型排列，而无序每次都再次对每个结构作用所有操作
                redundant = set()

                # loop over configurations
                for atoms_mark in atoms_gen(arg_sites):
                    arr_atoms_mark = numpy.array(atoms_mark)
                    flag_super = self._is_super(arr_atoms_mark, trans)
                    ahash = hash_atoms(atoms_mark)
                    if (ahash in redundant) or flag_super:
                        continue
                    else:
                        for p in perms:
                            atoms_transmuted = arr_atoms_mark[p]
                            redundant.add(hash_atoms(atoms_transmuted))

                        atoms = self._mark_to_atoms(arr_atoms_mark, sites)
                        # print("{:}".format(h.flatten()) +
                        #       '  ' + str(atoms))

                        yield Cell(supercell.lattice, supercell.positions, atoms_mark)

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

    def _is_super(self, arr_atoms_mark, translations):
        for t in translations[1:]:
            transmuted = arr_atoms_mark[t]
            if numpy.array_equal(transmuted, arr_atoms_mark):
                return True
        return False

    def _mark_to_atoms(self, arr_mark, sites):
        num_of_site_groups = len(sites)
        arr_atoms = arr_mark.reshape(num_of_site_groups, -1)
        # import pdb; pdb.set_trace()
        atoms = numpy.zeros_like(arr_atoms)
        for i, row in enumerate(arr_atoms):
            for j, v in enumerate(row):
                atoms[i][j] = sites[i][v]

        return atoms.flatten().tolist()

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
            hfpg = HFPG(self._pcell, h)
            # 此处的平移操作不必每次重新计算，因为相同snf平移操作相同
            # 可以用字典来标记查询。若没有这样的操作，那么就没有snf带来的效率提升。
            # For hnf with same snf, translations are same.
            quotient = hfpg.get_quotient()
            if not quotient in dict_trans:
                dict_trans[quotient] = hfpg.get_pure_translations(symprec)
            trans = dict_trans[quotient]

            supercell = self._pcell.extend(h)

            # 产生所有可能操作的置换操作
            # perm = hfpg.get_symmetry()
            rots = hfpg.get_pure_rotations(symprec)
            perms = self._get_perms_from_rots_and_trans(rots, trans)

            # TODO: 加入一个机制，来清晰的设定位点上无序的状态
            arg_sites_pcell = [len(i) for i in sites]
            arg_sites = numpy.repeat(arg_sites_pcell, volume)

            # redundant configurations do not want see again
            # 用于记录在操作作用后已经存在的构型排列，而无序每次都再次对每个结构作用所有操作
            redundant = set()

            # deg_total = 0
            # loop over configurations
            for atoms_mark in atoms_gen(arg_sites):
                arr_atoms_mark = numpy.array(atoms_mark)
                ahash = hash_atoms(atoms_mark)
                if (ahash in redundant):
                    continue
                else:
                    list_all_transmuted = []
                    for p in perms:
                        atoms_transmuted = arr_atoms_mark[p]
                        redundant.add(hash_atoms(atoms_transmuted))
                        # degeneracy
                        list_all_transmuted.append(atoms_transmuted)

                    arr_all_transmuted = numpy.array(list_all_transmuted)
                    deg = numpy.unique(arr_all_transmuted, axis=0).shape[0]

                    atoms = self._mark_to_atoms(arr_atoms_mark, sites)
                    # print("{:}".format(h.flatten()) +
                    #       '  ' + str(atoms) + '  ' + str(deg))
                    c = Cell(supercell.lattice, supercell.positions, atoms)
                    yield (c, deg)

            #     deg_total += deg
            # print(deg_total)

    def cons_specific_cell(self, sites, symprec=1e-5):
        """
        """
        lat_cell = self._cell.lattice
        lat_pcell = self._pcell.lattice
        # import pdb; pdb.set_trace()
        # TODO: use is_int_np_array
        mat = numpy.matmul(lat_cell, numpy.linalg.inv(lat_pcell))
        mat = mat.astype('intc')
        hfpg = HFPG(self._pcell, mat)

        # 产生所有可能操作的置换操作
        # perm = hfpg.get_symmetry()
        perms = hfpg.get_symmetry(symprec)

        # # TODO: 加入一个机制，来清晰的设定位点上无序的状态
        # arg_sites = [len(i) for i in sites]

        for c, d in self._remove_redundant(sites, perms):
            yield (c, d)

    def _remove_redundant(self, sites, perms):
        # TODO: 加入一个机制，来清晰的设定位点上无序的状态
        arg_sites = [len(i) for i in sites]
        # redundant configurations do not want see again
        # 用于记录在操作作用后已经存在的构型排列，而无序每次都再次对每个结构作用所有操作
        redundant = set()

        # deg_total = 0
        # loop over configurations
        for atoms_mark in atoms_gen(arg_sites):
            arr_atoms_mark = numpy.array(atoms_mark)
            ahash = hash_atoms(atoms_mark)
            if (ahash in redundant):
                continue
            else:
                list_all_transmuted = []
                for p in perms:
                    atoms_transmuted = arr_atoms_mark[p]
                    redundant.add(hash_atoms(atoms_transmuted))
                    # degeneracy
                    list_all_transmuted.append(atoms_transmuted)

                arr_all_transmuted = numpy.array(list_all_transmuted)
                deg = numpy.unique(arr_all_transmuted, axis=0).shape[0]

                atoms = self._mark_to_atoms(arr_atoms_mark, sites)
                # print("{:}".format(mat.flatten()) +
                #       '  ' +str(atoms) + '  ' + str(deg))
                # print(str(atoms) + '  ' + str(deg))

                c = Cell(self._cell.lattice, self._cell.positions, atoms)
                yield (c, deg)

        #     deg_total += deg
        # print(deg_total)
