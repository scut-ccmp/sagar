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
def configurations_nonredundant(pcell, sites, max_volume=2, symprec=1e-5):
    """
    parameters:

    pcell: Cell object, The primitive cell to be extended
    sites: disorderd sites infomation.

    yield:

    Cell object, a list of non-redundant configurations

    """
    # 该函数产生所有构型用于确定基态相图
    for v in range(1, max_volume + 1):
        hnfs = non_dup_hnfs(pcell, v, symprec)
        for h in hnfs:
            hfpg = HFPG(pcell, h)
            # TODO: 此处的平移操作不必每次重新计算，因为相同snf平移操作相同
            # 可以用字典来标记查询。若没有这样的操作，那么就没有snf带来的效率提升。
            translations = hfpg.get_pure_translations(symprec)
            supercell = pcell.extend(h)
            # TODO: 加入一个机制，来清晰的设定位点上无序的状态
            perm = hfpg.get_symmetry()
            arg_sites_pcell = [len(i) for i in sites]
            arg_sites = numpy.repeat(arg_sites_pcell, v)

            # actual results
            res = []
            # redundant configurations do not want see again
            # 用于记录在操作作用后已经存在的构型排列，而无序每次都再次对每个结构作用所有操作
            redundant = set()

            # loop over configurations
            for atoms_mark in atoms_gen(arg_sites):
                # import pdb; pdb.set_trace()
                arr_atoms_mark = numpy.array(atoms_mark)
                # import pdb; pdb.set_trace()
                flag_super = _is_super(arr_atoms_mark, translations)
                ahash = hash_atoms(atoms_mark)
                # import pdb; pdb.set_trace()
                if (ahash in redundant) or flag_super:
                    continue
                else:
                    for p in perm:
                        atoms_transmuted = arr_atoms_mark[p]
                        redundant.add(hash_atoms(atoms_transmuted))

                    atoms = _mark_to_atoms(arr_atoms_mark, sites)
                    print("{:}".format(h.flatten()) +
                          '  ' + str(atoms))

                    # c = Cell(supercell.lattice, supercell.positions, atoms)
                    # yield Cell(supercell.lattice, supercell.positions, atoms_mark)
                    # print(c)
                    # yield c


def _is_super(arr_atoms_mark, translations):
    for t in translations[1:]:
        transmuted = arr_atoms_mark[t]
        if numpy.array_equal(transmuted, arr_atoms_mark):
            return True
    return False

def _mark_to_atoms(arr_mark, sites):
    num_of_site_groups = len(sites)
    arr_atoms = arr_mark.reshape(num_of_site_groups, -1)
    # import pdb; pdb.set_trace()
    atoms = numpy.zeros_like(arr_atoms)
    for i, row in enumerate(arr_atoms):
        for j, v in enumerate(row):
            atoms[i][j] = sites[i][v]

    return atoms.flatten().tolist()
