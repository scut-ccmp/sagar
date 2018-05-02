# coding: utf-8
# Distributed under the terms of the MIT license.
import os
import json
import xxhash
# from hat_trie import Trie

import numpy as np
import spglib

from ababe.stru.element import Specie

# Loads bucky-ball structure data from json file
with open(os.path.join(os.path.dirname(__file__),
                            "buckyball.json"), "rt") as f:
    _buckyball = json.load(f)

# Loads buckyball's permutation table from .npy
_perm_table = np.load(os.path.join(os.path.dirname(__file__),
                                    "bucky_permu_table.npy"))

class Structure(object):
    """
    Class to generate bucky-ball related structure parameters.
    """

    def __init__(self, numbers):
        self._lattice = np.array(_buckyball["lattice"])
        numbers = np.array(numbers)

        # Sorting positions (x,y,z)
        init_positions = np.array(_buckyball["positions"])
        init_index = self._get_new_id_seq(init_positions, numbers)
        self._positions = init_positions[init_index]

        self._atom_numbers = numbers
        self._spg_cell = (self._lattice, self._positions, self._atom_numbers)
        self._carbon = Specie("C")
        self._sym_perm = list(_perm_table)

    @property
    def spg_cell(self):
        return self._spg_cell

    def get_spacegroup(self):
        return spglib.get_spacegroup(self._spg_cell, symprec=1e-4)

    def get_speckle_num(self, sp):
        from collections import Counter
        atom = sp.Z
        num_count = Counter(self._atom_numbers)
        num = num_count[atom]
        #num = self._atom_numbers.count(atom)
        return num

    def get_symmetry(self):
        """
        Symmetry operations are obtained as a dictionary. 
        The key rotation contains a numpy array of integer, 
        which is “number of symmetry operations” x “3x3 matrices”. 
        The key translation contains a numpy array of float, 
        which is “number of symmetry operations” x “vectors”. 
        """
        symmetry = spglib.get_symmetry(self._spg_cell, symprec=1e-4)
        return symmetry

    def get_symmetry_permutation(self):
        """
        This a object function to get the permutation group operators.
        Represented as a table.
        """
        # sym_perm = []
        # numbers = [i for i in range(60)]
        # sym_mat = spglib.get_symmetry(self._spg_cell, symprec=1e-4)
        # ops = [(r,t) for r, t in zip(sym_mat['rotations'], sym_mat['translations'])]
        # for r, t in ops:
        #     pos_new = np.transpose(np.matmul(r, np.transpose(self._positions))) + t
        #     perm = self._get_new_id_seq(pos_new, numbers)
        #     sym_perm.append(perm)

        return self._sym_perm

    @staticmethod
    def _get_new_id_seq(pos, numbers):
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

    def get_name(self):
        """
        For the reason all structure have same lattice and positions
        only atom sequences are diff. Therefore as the hashable name
        and dict key.
        """
        return self._atom_numbers.tostring()

    def to_gen(self):
        """
        This a method convert a number seqents to a one element
        generator.
        """
        l = [self._atom_numbers]
        g = (n for n in l)
        return g

    @staticmethod
    def gen_speckle(sp, noa):
        """
        This method creates speckle structures which have speckles 
        number = noa.
        """
        i_sea = 6
        i_speckle = sp.Z
        from itertools import combinations
        for comb_index in combinations(range(60), noa):
            atom_numbers = [i_sea]*60
            for index in comb_index:
                atom_numbers[index] = i_speckle
            yield np.array(atom_numbers)

    @staticmethod
    def help_add_one_speckle(l, sp):
        atom = sp.Z
        for index, val in enumerate(l):
            l_new= list(l)
            if val != atom:
                l_new[index] = atom
                yield l_new

    def add_one_speckle_generator(self, gen, sp):
        """
        input a structure generator(mostly nonduplicate)
        output a generator with one more speckle.
        This a method give duplicate structures which have one more speckle
        than the input structures.
        """
        atom = sp.Z
        idy_seq = dict()
        for s_atoms in gen:
            for index, val in enumerate(s_atoms):
                arr_new = s_atoms.copy()
                #print(arr_new)
                if val != atom:
                    arr_new[index] = atom
                    arr_idy = self._get_atom_seq_identifier(arr_new)
                    #arr_idy = xxhash.xxh32(arr_new).intdigest()
                    if arr_idy not in idy_seq:
                        yield arr_new
                        idy_seq[arr_idy] = None

    @staticmethod
    def _get_atom_seq_identifier(numbers):
        """
        This method convert a list to a immutable string, which used
        as an identifier of diffrent structures, can be move to
        outerside class.
        """
        num_hash = xxhash.xxh32(numbers).intdigest()
        return num_hash

    def _update_isoset(self, isoset, atoms, sym_perm):
        for ind in sym_perm:
            atoms_new = atoms[ind]
            id_stru = self._get_atom_seq_identifier(atoms_new)
            isoset[id_stru] = None

        #return isoset

    def to_nodup_generator(self, dup_gen):
        """
        input: a generator with duplicate structures
        output: a generator with no structures dupicated
        This a method filter the duplicate structure to nonduplicates.
        """
        sym_perm = self.get_symmetry_permutation()

        isoset = dict()
        for atoms in dup_gen:
            id_stru = self._get_atom_seq_identifier(atoms)
            if id_stru not in isoset:
                yield atoms
                self._update_isoset(isoset, atoms, sym_perm)

    def nodup_generator(self):

        while True:
            yield atoms

    @staticmethod
    def all_speckle_gen(bucky_stru, n_max, sp):
        from itertools import tee
        gen = bucky_stru.to_gen()
        out_gen, gen = tee(gen, 2)
        yield out_gen
        n_init = bucky_stru.get_speckle_num(sp)
        for i in range(n_init, n_max):
            gen = bucky_stru.add_one_speckle_generator(gen, sp)
            gen = bucky_stru.to_nodup_generator(gen)

            out_gen, gen = tee(gen, 2)
            yield out_gen

