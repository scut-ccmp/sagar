# coding: utf-8
# Distributed under the terms of the MIT License.

import pdb
import sys
import os.path

from ababe.stru.scaffold import SitesGrid, CStru, GeneralCell
from ababe.stru.element import GhostSpecie, Specie
from itertools import combinations, zip_longest
from progressbar import ProgressBar

import numpy as np
import spglib
from spglib import get_symmetry
import os
import xxhash
from itertools import tee
# Filename sogen is for Site-Occupy-GENerator


def numbers2id(numbers_arr):
    num_hash = xxhash.xxh64(numbers_arr).intdigest()
    return num_hash


class OccupyGenerator(object):
    """
    The class takes a GeneralCell instance as input,
    and can produces a generator GeneralCell instances which are
    nonduplicated to each others.
    """
    def __init__(self, general_cell):
        self.init_cell = general_cell
        self.lattice = general_cell.lattice
        self.positions = general_cell.positions
        self.numbers = general_cell.numbers

        self.symmetry_permutation = general_cell.get_symmetry_permutation()
        self.wyckoffs = general_cell.get_wyckoffs()

    def is_equivalent(self, cell_i, cell_other):
        numbers_i = cell_i.numbers
        numbers_other = cell_other.numbers
        for perm in self.symmetry_permutation:
            new_numbers = numbers_i[perm]
            if np.array_equal(new_numbers, numbers_other):
                return True

        return False

    def gen_dup_unitary(self, n, sp):
        init_numbers = self.init_cell.numbers
        num_count = init_numbers.size   # number of atoms in structure
        i_speckle = sp.Z
        from itertools import combinations
        for comb_index in combinations(range(num_count), n):
            numbers = init_numbers.copy()
            for index in comb_index:
                numbers[index] = i_speckle
            yield GeneralCell(self.lattice, self.positions, numbers)

    def gen_nodup_unitary(self, n, sp):
        dup = self.gen_dup_unitary(n, sp)
        # sym_perm = self.symmetry_permutation

        # isoset = dict()
        # for cell in dup:
        #     cell_id = cell.id
        #     if cell_id not in isoset:
        #         yield cell
        #         # from ababe.stru.io import VaspPOSCAR
        #         # pdb.set_trace()
        #         self._update_isoset(isoset, cell.numbers, sym_perm)
        return self.gen_2nodup_gen(dup)

    def gen_dup(self, wy, n, sp):
        init_numbers = self.init_cell.numbers
        i_speckle = sp.Z
        wyckoffs = self.wyckoffs
        sp_ind = [i for i, w in enumerate(wyckoffs) if w is wy]
        # pdb.set_trace()
        for comb_index in combinations(sp_ind, n):
            numbers = init_numbers.copy()
            for index in comb_index:
                numbers[index] = i_speckle
            yield GeneralCell(self.lattice, self.positions, numbers)

    def gen_nodup(self, wy, n, sp):
        dup = self.gen_dup(wy, n, sp)
        return self.gen_2nodup_gen(dup)

    def gen_dup_of_ele(self, ele, n, sp):
        # ele is a Specie instance

        ele_num = ele.Z
        init_numbers = self.init_cell.numbers
        i_speckle = sp.Z
        sp_ind = [i for i, e in enumerate(init_numbers) if e == ele_num]

        for comb_index in combinations(sp_ind, n):
            numbers = init_numbers.copy()
            for index in comb_index:
                numbers[index] = i_speckle
            yield GeneralCell(self.lattice, self.positions, numbers)

    def gen_nodup_of_ele(self, ele, n, sp):
        dup = self.gen_dup_of_ele(ele, n, sp)
        return self.gen_2nodup_gen(dup)

    def gen_dup_exch(self, sp1, sp2, n):
        init_numbers = self.init_cell.numbers
        n_sp1 = sp1.Z
        n_sp2 = sp2.Z
        sp1_ind = [i for i, e in enumerate(init_numbers) if e == n_sp1]
        sp2_ind = [i for i, e in enumerate(init_numbers) if e == n_sp2]

        for ex_sp1 in combinations(sp1_ind, n):
            for ex_sp2 in combinations(sp2_ind, n):
                numbers = init_numbers.copy()
                for ind_sp1, ind_sp2 in zip(ex_sp1, ex_sp2):
                    numbers[ind_sp1] = n_sp2
                    numbers[ind_sp2] = n_sp1
                yield GeneralCell(self.lattice, self.positions, numbers)


    def gen_nodup_exch(self, sp1, sp2, n):
        dup = self.gen_dup_exch(sp1, sp2, n)
        return self.gen_2nodup_gen(dup)

    def _gen_dup_trinary_alloy(self, sp1, n1, sp2, n2):
        init_numbers = self.init_cell.numbers
        isp1 = sp1.Z
        isp2 = sp2.Z

        sp_ind_origin = [i for i, s in enumerate(init_numbers)]
        for sp1_comb_index in combinations(sp_ind_origin, n1):
            sp_ind_bin = [x for x in sp_ind_origin if x not in sp1_comb_index]
            for sp2_comb_index in combinations(sp_ind_bin, n2):
                numbers = init_numbers.copy()
                for i1, i2 in zip_longest(sp1_comb_index, sp2_comb_index):
                    if i1 is not None:
                        numbers[i1] = isp1
                    if i2 is not None:
                        numbers[i2] = isp2
                yield GeneralCell(self.lattice, self.positions, numbers)
                # pdb.set_trace()

    def gen_nodup_trinary_alloy(self, sp1, n1, sp2, n2):
        dup = self._gen_dup_trinary_alloy(sp1, n1, sp2, n2)
        return self.gen_2nodup_gen(dup)

    def gen_add_one_speckle_unitary(self, gen, sp):
        """
        input a structure generator __ spg_cell(mostly nonduplicate)
        output a generator with one more speckle.
        This a method give duplicate structures which have one more speckle
        than the input structures.
        """
        atom = sp.Z
        id_db = dict()
        for cell in gen:
            for index, val in enumerate(cell.numbers):
                numbers_new = cell.numbers.copy()
                if atom != val:
                    numbers_new[index] = atom
                    num_id = numbers2id(numbers_new)
                    if num_id not in id_db:
                        yield GeneralCell(cell.lattice, cell.positions,
                                          numbers_new)
                        id_db[num_id] = None

    def gen_add_one_speckle(self, gen, wy, sp):
        atom = sp.Z
        id_db = dict()
        for cell in gen:
            wy_ele = zip(self.wyckoffs, cell.numbers)
            for index, tuple_w_e in enumerate(wy_ele):
                numbers_new = cell.numbers.copy()
                if tuple_w_e[0] == wy and atom != tuple_w_e[1]:
                    numbers_new[index] = atom
                    num_id = numbers2id(numbers_new)
                    if num_id not in id_db:
                        yield GeneralCell(cell.lattice, cell.positions,
                                          numbers_new)
                        id_db[num_id] = None

    def gen_add_one_speckle_of_ele(self, gen, ele, sp):
        ele_Z = ele.Z
        atom = sp.Z
        id_db = dict()
        for cell in gen:
            # site_ele = zip(self.wyckoffs, cell.numbers)
            for index, e in enumerate(cell.numbers):
                numbers_new = cell.numbers.copy()
                if e == ele_Z and e != atom:
                    numbers_new[index] = atom
                    num_id = numbers2id(numbers_new)
                    if num_id not in id_db:
                        yield GeneralCell(cell.lattice, cell.positions,
                                          numbers_new)
                        id_db[num_id] = None

    def gen_2nodup_gen(self, dup_gen):
        sym_perm = self.symmetry_permutation

        isoset = dict()
        bar = ProgressBar()
        for cell in bar(dup_gen):
            if cell.id not in isoset:
                ## TODO: combine degeneracy here
                yield cell
                self._update_isoset(isoset, cell.numbers, sym_perm)

    @staticmethod
    def _update_isoset(isoset, numbers, sym_perm):
        for ind in sym_perm:
            numbers_new = numbers[ind]
            cell_id = numbers2id(numbers_new)
            isoset[cell_id] = None

    def all_speckle_gen_unitary(self, n, sp):
        gen = (i for i in [self.init_cell])
        n_init = self.init_cell.get_speckle_num(sp)
        print("Mission: Replace with {0:4}, up to\
               {1:4d}...".format(sp.name, n))
        for i in range(n_init, n):
            gen = self.gen_add_one_speckle_unitary(gen, sp)
            gen = self.gen_2nodup_gen(gen)

            out_gen, gen = tee(gen, 2)
            yield out_gen

    def all_speckle_gen(self, n, wy, sp):
        gen = (i for i in [self.init_cell])
        n_init = self.init_cell.get_speckle_num(sp)
        print("Mission: Replace with {0:4}, up to\
              {1:4d}, in wyckoff site {2:3}...".format(sp.name, n, wy))
        for i in range(n_init, n):
            gen = self.gen_add_one_speckle(gen, wy, sp)
            gen = self.gen_2nodup_gen(gen)

            out_gen, gen = tee(gen, 2)
            yield out_gen

    def all_speckle_gen_of_ele(self, n, ele, sp):
        gen = (i for i in [self.init_cell])
        n_init = self.init_cell.get_speckle_num(sp)
        print("Mission: Replace with {0:4}, up to\
              {1:4d}, in wyckoff site {2:3}...".format(sp.name, n, ele.Z))
        for i in range(n_init, n):
            gen = self.gen_add_one_speckle_of_ele(gen, ele, sp)
            gen = self.gen_2nodup_gen(gen)

            out_gen, gen = tee(gen, 2)
            yield out_gen


def is_stru_equal(struA, struB, ops):
    bA, posA, atom_numA = struA.get_cell()
    bB, posB, atom_numB = struB.get_cell()
    id_struA = _get_id_seq(posA, atom_numA)

    is_equal = False
    for r, t in ops:
        pos_new = np.transpose(np.matmul(r, np.transpose(posB))) + t
        id_stru = _get_id_seq(pos_new, atom_numB)
        if id_stru == id_struA:
            is_equal = True

    return is_equal

def _get_id_seq(pos, arr_num):

    # from fractions import Fraction
    # transfer the atom position into >=0 and <=1
    pos = np.around(pos, decimals=10)
    func_tofrac = np.vectorize(lambda x: round((x % 1), 3))
    o_pos = func_tofrac(pos)
    # round_o_pos = np.around(o_pos, decimals=3)
    # z, y, x = round_o_pos[:, 2], round_o_pos[:, 1], round_o_pos[:, 0]
    z, y, x = o_pos[:, 2], o_pos[:, 1], o_pos[:, 0]
    ind_sort = np.lexsort((z, y, x))
    id_seq = str(arr_num[ind_sort])

    return id_seq

def _get_atom_seq_identifier(numbers):
    return str(list(numbers))

def _update_isoset(isoset, numbers, sym_perm):
    for ind in sym_perm:
        # pos_new = np.transpose(np.matmul(r, np.transpose(pos))) + t
        sequence_new = numbers[ind]
        id_stru = _get_atom_seq_identifier(sequence_new)
        # isoset_cstru.add(id_stru)
        # isoset.update(isoset_cstru)
        isoset.add(id_stru)

    return isoset

def get_new_id_seq(pos, numbers):
    """
    A helper function to produce the new sequence of the transformed
    structure. Algs is sort the position back to init and use the index
    to sort numbers.
    """
    # transfer the atom position into >=0 and <=1
    pos = np.around(pos, decimals=5)
    func_tofrac = np.vectorize(lambda x: round((x % 1), 3))
    o_pos = func_tofrac(pos)
    # round_o_pos = np.around(o_pos, decimals=3)
    # z, y, x = round_o_pos[:, 2], round_o_pos[:, 1], round_o_pos[:, 0]
    z, y, x = o_pos[:, 2], o_pos[:, 1], o_pos[:, 0]
    inds = np.lexsort((z, y, x))

    return inds

def get_permutation_cell(cell):
    lat, pos, num = cell
    atom_num = len(cell[2])
    numbers = [i for i in range(atom_num)]
    sym = get_symmetry(cell, symprec=1e-3)
    ops = [(r, t) for r, t in zip(sym['rotations'], sym['translations'])]
    sym_perm = []
    for r,t in ops:
        pos_new = np.transpose(np.matmul(r, np.transpose(pos))) + t
        perm = get_new_id_seq(pos_new, numbers)
        sym_perm.append(perm)

    return sym_perm

def gen_nodup_cstru(lattice, sea_ele, size, speckle, num):
    d, w, l = size
    ele_sea = SitesGrid.sea(d, w, l, sea_ele)
    cell_mother_stru = CStru(lattice, ele_sea).get_cell()

    sym_perm = get_permutation_cell(cell_mother_stru)

    # For testing: Show that the first unit matrix convert to range(num) perm_operator
    # print(sym_perm[0])

    gen_dup_cstrus = CStru.gen_speckle(lattice, sea_ele, size, speckle, num)

    # Add the progress bar when looping
    from scipy.special import comb
    number_of_structures = comb((d*w*l), num)
    bar = ProgressBar(max_value=number_of_structures)

    isoset = set()
    for cstru in bar(gen_dup_cstrus):
        b, pos, atom_num = cstru.get_cell()
        #id_cstru = _get_id_seq(pos, atom_num)
        id_cstru = _get_atom_seq_identifier(atom_num)
        # print(id_cstru)
        if id_cstru not in isoset:
            # print(len(sym_perm))
            # print(len(isoset))
            # print(cstru.get_array())
            yield cstru
            _update_isoset(isoset, atom_num, sym_perm)

def default(str):
    return str + '  [Default: %(default)s]'

def lat_dict(lattice):
    from math import sqrt
    lat = { 'bcc': [[-0.5, -0.5, -0.5],
                [-0.5,  0.5,  0.5],
                [ 0.5, -0.5,  0.5]],
            'fcc': [[0, 0.5, 0.5],
                    [0.5, 0, 0.5],
                    [0.5, 0.5, 0]],
            'scc': [[1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1]],
            'triflat': [[0, 0, 20],
                        [1, 0, 0],
                        [0.5, sqrt(3)/2, 0]]
            }

    return lat[lattice]

# This function is used for remove the structures conflict with
# the defined restricted condition
# input: a generator to produce structures
# output: a generator of structures satisfied with the restricted
# condition.
def is_speckle_disjunct(cstru, speckle):
    m = cstru.m
    sites_arr = cstru.get_array()
    ele = speckle.Z

    pool_sites_arr = _pool_sites(sites_arr)
    ele_index = np.argwhere(pool_sites_arr==ele)
    points = np.array([_index2coor(ind, m) for ind in ele_index])
    min_d = _min_dist(points)
    is_disjunct = min_d > 1.01

    return is_disjunct

def _min_dist(points):
    # get the closest pair of points
    # Brute-force algorithm
    min_distance = 9999
    pairs = combinations(points, 2)
    for pair in pairs:
        if _dist(pair[0], pair[1]) < min_distance:
            min_distance = _dist(pair[0], pair[1])

    return min_distance

def _dist(p,q):
    dist = np.linalg.norm(p-q)
    return dist

def _pool_sites(sites_arr):
    d, w, l = np.shape(sites_arr)
    pool = sites_arr
    # pool the elements of outer dimension (depth)
    depth_d = pool[0, :, :].reshape(1,w,l)
    pool = np.concatenate((pool, depth_d), axis=0)
    # pool the elements of meddle dimension (width)
    width_d = pool[:, 0, :].reshape(d+1, 1, l)
    pool = np.concatenate((pool, width_d), axis=1)
    # pool the elements of inner dimension (length)
    length_d = pool[:, :, 0].reshape(d+1, w+1, 1)
    pool = np.concatenate((pool, length_d), axis=2)

    return pool

def _index2coor(ind, m):
    m_arr = np.array(m)
    d, w, l = ind
    v1 = m_arr[0]*d
    v2 = m_arr[1]*w
    v3 = m_arr[2]*l
    cood = np.array([v1[0]+v2[0]+v3[0], v1[1]+v2[1]+v3[1], v1[2]+v2[2]+v3[2]])
    return cood

# def main():
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--lattice', choices=['bcc', 'fcc', 'scc', 'triflat'],
    #                         help='Lattice type of grid conventional cell')
    # parser.add_argument('-b', '--base', dest='sea', required=True,
    #                         help='Element abbreviation of the base specie')
    # parser.add_argument('-g', '--size', nargs=3, dest='size', required=True,
    #                         help='Grid size of structure', type=int)
    # parser.add_argument('-s', '--speckle', dest='speckle', required=True,
    #                         help='Element abbreviation of the speckle specie')
    # parser.add_argument('-n', '--num', dest='number', type=int,
    #                         help=default('Number of speckles filled in the base'), default=2)
    # parser.add_argument('-o', '--output-type',
    #                         help=default('Output type of generated non-duplicate periodic grid structure'),
    #                             default='normal')

    # args = parser.parse_args()
    # size = args.size
    # sea_ele = Specie(args.sea)
    # speckle = Specie(args.speckle)

    # nodup_gen = gen_nodup_cstru(lat_dict(args.lattice), sea_ele, size, speckle, args.number)
    # with open('allstru.txt', 'w') as f:
    #     for s in nodup_gen:
    #         basis, pos, atom = s.get_cell()
    #         f.write('ONE NEW STRUCTURE:\n')
    #         f.write('The basis is:\n')
    #         f.write('\n'.join(str(line) for line in basis))

    #         f.write('\nThe position is:\n')
    #         f.write('\n'.join(str(line) for line in pos))

    #         f.write('\nThe elements is:\n')
    #         f.write(str(atom))

    #         f.write('\n\n\n')
