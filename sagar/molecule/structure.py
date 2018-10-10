# -*- coding: utf-8 -*-
"""
Created on Sat Aug 25 16:41:45 2018
@author: hecc

Modified by: Jason Yu
"""
import numpy

from sagar.molecule import symmetry
from sagar.toolkit.mathtool import closest_pair
from sagar.element.base import periodic_table_dict, get_symbol, symbol2number


class Molecule(object):
    """
    Molecule object represent a 1D Molecule structure.

    parameters:

    positions: n tuples(x,y,z) in fraction. 可多态初始化
    atoms: list of atoms, can be atomic number (int), can be
    atomic symbol (string),
    represent atom in periodic table. (用全称初始化？？谁TM会这么用？)
    """

    def __init__(self,  positions, atoms):
        # TODO: magmoms init setting
        self._atom_numbers = len(atoms)
        self._positions = numpy.array(positions).reshape((-1, 3))
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

    @property
    def positions(self):
        return self._positions

    @property
    def atoms(self):
        return self._atoms

    def __repr__(self):
        def _repr(number):
            return "{:9.6f}".format(number)

        sites = zip(self._positions, self._atoms)
        out_pos = []
        out_pos.append("Sites:")
        for s in sites:
            o = ' '.join(map(_repr, s[0])) + ' ' + get_symbol(s[1])
            out_pos.append(o)
        outs = out_pos
        return "\n".join(outs)

    def get_symmetry_permutation(self, symprec=1e-3):

        return symmetry.get_permutations(self._positions,
                                         self._atoms.tolist(),
                                         symprec)

    def check(self, elements=None, limit=0.1, warn=False):
        """
        该方法用于自查对象中的位点是否过近
        若过近则抛出一个warning
        """
        if elements is None:
            mol = self
        else:
            # 选取要check的元素，构建新的胞
            ele_num = [symbol2number(s) for s in elements]
            positions = []
            atoms = []
            for idx, e in enumerate(self.atoms):
                if e in ele_num:
                    positions.append(self.positions[idx])
                    atoms.append(self.atoms[idx])
            mol = self.__class__(positions, atoms)

        points = mol.positions
        if closest_pair(points) < limit:
            if warn is True:
                import warnings
                warnings.warn("some atoms are too close(< {:f}), check cell".format(
                    limit), RuntimeWarning)
            return False
        else:
            return True
