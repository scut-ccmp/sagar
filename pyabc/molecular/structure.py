#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 25 16:41:45 2018

@author: hecc
"""

import numpy
import Molecular_symmetry


"""
periodic_table_dict
Thanks for Yanan Wu. Email:
"""
periodic_table_dict = {'Vacc': 0,
                       'H': 1, 'He': 2,
                       'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
                       'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
                       'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
                       'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54,
                       'Cs': 55, 'Ba': 56, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86,
                       'Fr': 87, 'Ra': 88}


def get_symbol(atom):
    """
    get_symbol return symbol of atomic number

    parameter:

    atom: int, atomic number.
    """
    for key, value in periodic_table_dict.items():
        if atom == value:
            return str(key)
    return "NaN_x"


def symbol2number(symbol):
    return periodic_table_dict[symbol]


class Cell(object):
    """
    Cell object represent a crystal structure.

    parameters:

    lattice: 3x3 1D-list, lattice of cell
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

    def get_symmetry_permutation(self, pres=1e-3):

        return Molecular_symmetry.get_permutations(self._positions, \
                                                   self._atoms.tolist(), pres)



    def check(self, limit=0.1):
        """
        该方法用于自查对象中的位点是否过近
        若过近则抛出一个warning
        """
        n = len(self.atoms)
        warn_list = []
        for ii in range(n):
            for jj in range(ii+1, n):
                temp = numpy.linalg.norm(self._positions[ii]-self._positions[jj])
                if temp < limit:
                    warn_list.append(str(ii)+'-'+str(jj))
        if warn_list:
            import warnings
            warnings.warn( "["+(', '.join(['%s']*len(warn_list))+"]")%tuple(\
                             warn_list)+
                          (" atom couples are too close { < %.2f } \n")%(limit,))
            return True

        else:
            return False
