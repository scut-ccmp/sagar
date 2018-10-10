# -*- coding: utf-8 -*-
import numpy
from itertools import product

from sagar.toolkit.mathtool import binomialCoeff

def remove_redundant(mol_positions, sites, perms, e_num=None, method='jshash'):
    if method == 'jshash':
        for i in remove_redundant_by_hash(mol_positions, sites, perms, e_num):
            yield i
    if method == 'ccsort':
        pass

def remove_redundant_by_hash(mol_positions, sites, perms, e_num):
    """
    输入一个分子坐标`mol_positions`
    和每个坐标位点上的可能取代情况`sites`
    母体结构所有的置换操作`perms`

    options:
        位点上元素的比例`e_num`
    """
    # perms = self._get_perms_from_rots_and_trans(rots, trans)
    # TODO: 加入一个机制，来清晰的设定位点上无序的状态
    sites = numpy.array(sites)
    arg_sites = [len(i) for i in sites]
    # redundant configurations do not want see again
    # 用于记录在操作作用后已经存在的构型排列，而无序每次都再次对每个结构作用所有操作
    redundant = set()

    # deg_total = 0
    # loop over configurations
    for atoms_mark in _atoms_gen(arg_sites, e_num):
        arr_atoms_mark = numpy.array(atoms_mark)
        ahash = _hash_atoms(atoms_mark)

        if ahash in redundant:
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

            atoms = _mark_to_atoms(arr_atoms_mark, sites)

            # 返回一个包含不重复位置和原子类别的tuple
            # TODO: 合并成为一个独特的类？？
            m = (mol_positions, atoms)
            yield (m, deg)

def _mark_to_atoms(arr_mark, sites):
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
        # import pdb; pdb.set_trace()
        disorder_site = [s for s in args if s > 1]
        num_disorder_site = len(disorder_site)
        if num_disorder_site != sum(e_num):
            raise ValueError("concentration given error, wanted sum {:d}, got {:d}".format(
                num_disorder_site, sum(e_num)))
        arr_arrange = _serial_int_to_arrangement(e_num)

        # 若该位点只有一个元素，则在arr_arrange中加入指定列
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
            n_color = a
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

            open_slots -= n_color
    return arr_arrangement
