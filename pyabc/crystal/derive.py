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


def configurations_nonredundant(pcell, sites, volume=1, symprec=1e-5):
    """
    parameters:

    pcell: Cell object, The primitive cell to be extended
    sites: disorderd sites infomation.

    yield:

    Cell object, a list of non-redundant configurations

    在较大体积的构型中，可能有结构已经在小的体积中出现过。

    该函数中没有comprec参数，是因为构型(configurations)之间的对比
    通过原子排列的对比实现，没有坐标的比较。
    """
    hnfs = non_dup_hnfs(pcell, volume, symprec)
    for h in hnfs:
        hfpg = HFPG(pcell, h)
        # TODO: 此处的平移操作不必每次重新计算，因为相同snf平移操作相同
        # 可以用字典来标记查询。若没有这样的操作，那么就没有snf带来的效率提升。
        # translations = hfpg.get_pure_translations(symprec)
        supercell = pcell.extend(h)
        # TODO: 加入一个机制，来清晰的设定位点上无序的状态
        perm = hfpg.get_symmetry()
        arg_sites_pcell = [len(i) for i in sites]
        arg_sites = numpy.repeat(arg_sites_pcell, volume)

        # actual results
        res = []
        # redundant configurations do not want see again
        # 用于记录在操作作用后已经存在的构型排列，而无序每次都再次对每个结构作用所有操作
        redundant = set()

        # loop over configurations
        for atoms_mark in atoms_gen(arg_sites):
            # import pdb; pdb.set_trace()
            atoms_mark_arr = numpy.array(atoms_mark)
            ahash = hash_atoms(atoms_mark)
            if ahash in redundant:
                continue
            else:
                res.append(atoms_mark)

            # for diff hnf loop over rots and trans, aka transfromations
            for p in perm:
                atoms_transmuted = atoms_mark_arr[p]
                redundant.add(hash_atoms(atoms_transmuted))

        print("hnf is:")
        print(h)
        print(supercell.lattice)
        print(supercell.positions)
        for atoms in res:
            # yield Cell(supercell.lattice, supercell.positions, atoms)
            print(atoms)
        print('\n')
