# -*- coding: utf-8 -*-
import collections
import operator
import numpy

from sagar.molecule.structure import Molecule
from sagar.element.base import get_symbol


def read_vasp(filename='XYZ'):
    """
    Import filename with .xyz suffix

    parameter:
    filename: string, the filename

    return: Molecule object.
    """
    # TODO: read velocities, now not supported. or not needed?
    with open(filename, "r") as f:
        # _read_string return Molecule object.
        return _read_string(f.read())


def _read_string(data):
    """
    _read_string make io easy to be tested.

    parameter: string of xyz input

    return: Molecule object
    """
    lines = [l for l in data.split('\n') if l.rstrip()]

    total_atoms = int(lines[0])

    comment = lines[1]

    positions = []
    atoms = []
    for i in range(2, 2 + total_atoms):
        s = lines[i].split()
        atoms.append(s[0])
        vec = float(s[1]), float(s[2]), float(s[3])
        positions.append(vec)

    return Molecule(positions, atoms)

def write_vasp(mol, filename='XYZ', suffix='.xyz', long_format=True):
    """
    write xyz type into file.
    always write atoms sorted xyz file.

    parameters:

    mol: Molecule object, the Molecule that you wanna write into xyz.
    filename: string, filename of output file, default='XYZ'
    suffix: string, suffix of filename, default='.xyz'
    long_format: bool, if True format %.16f will be write, else %.6f
    ref: https://gitlab.com/ase/ase/blob/master/ase/io/vasp.py

    if optional parameters (filename and suffix) are not set,
    the filename will be 'XYZ.xyz'
    """
    filname_suffix = ''.join([filename, suffix])
    with open(filname_suffix, "w") as f:
        f.write(_write_string(cell, long_format))

def _write_string(mol, long_format, print_vacc=False):
    """
    _write_string make io easy to be tested.

    return: string represtent XYZ file fomat
    """

    # 对原子种类合并排序，用以产生体系名称和原子顺序数目和正确的坐标排序
    # sorted is a list of tuple(atom, na)
    atoms_dict = collections.Counter(mol.atoms)
    if not print_vacc:
        del atoms_dict[0]
    sorted_symbols = sorted(atoms_dict.items(), key=operator.itemgetter(0))

    list_symbols = ["{:}{:}".format(get_symbol(atom), na)
                    for atom, na in sorted_symbols]

    total_atoms = 0
    for n in atoms_dict.values():
        total_atoms += n
    total_atoms = str(total_atoms)
    total_atoms += '\n'

    comment = ' '.join(list_symbols)
    comment += '\n'

    # argsort atoms and resort coor
    idx = numpy.argsort(mol.atoms)
    coord = mol.positions[idx]
    atoms = mol.atoms[idx]
    positions_string = ""
    if long_format:
        pos_form = '19.16f'
    else:
        pos_form = '9.6f'

    for i, vec in enumerate(coord):
        if atoms[i] == 0:
            continue
        positions_string += ' ' + get_symbol(atoms[i])
        for v in vec:
            positions_string += '{:{form}}'.format(v, form=pos_form)
        positions_string += '\n'

    xyz_string = ''.join([total_atoms,
                          comment,
                          positions_string])
    return xyz_string
