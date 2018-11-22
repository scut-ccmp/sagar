# -*- coding: utf-8 -*-
import numpy
import warnings
import operator
import collections

from sagar.crystal.structure import Cell
from sagar.element.base import get_symbol


def read_vasp(filename='POSCAR'):
    """
    Import POSCAR/CONTCAR or filename with .vasp suffix

    parameter:
    filename: string, the filename

    return: Cell object.
    """
    # TODO: read velocities, now not supported. or not needed?
    with open(filename, "r") as f:
        # _read_string return Cell object.
        return _read_string(f.read())


def _read_string(data):
    """
    _read_string make io easy to be tested.

    parameter: string of vasp input

    return: Cell object
    """
    lines = [l for l in data.split('\n') if l.rstrip()]

    name = lines[0]

    lattice_scale = float(lines[1].split()[0])

    # lattice vectors
    lattice = []
    for i in [2, 3, 4]:
        s = lines[i].split()
        vec = float(s[0]), float(s[1]), float(s[2])
        lattice.append(vec)
    lattice = numpy.array(lattice)

    if lattice_scale < 0:
        # In vasp , a negative scale factor is treated as a volume.
        # http://pymatgen.org/_modules/pymatgen/io/vasp/inputs.html#POSCAR
        vol = abs(numpy.linalg.det(lattice))
        lattice *= (-lattice_scale / vol) ** (1 / 3)
    else:
        lattice *= lattice_scale

    # atoms
    vasp5 = False
    _fifth_line = lines[5].split()
    # VASP 5.x use the fifth line to represent atomic symbols
    try:
        for i in _fifth_line:
            int(i)
        numofatoms = _fifth_line
    except ValueError:
        vasp5 = True
        atomtypes = _fifth_line
        numofatoms = lines[6].split()  # list of string here

    if not vasp5:
        warnings.warn("symbols of elements in fifth line are missing, "
                      "all atoms are init to NaN_i (i=0,1,2...)", UserWarning, stacklevel=2)
        atomtypes = [str("NaN_{:}".format(i)) for i in range(len(numofatoms))]

    atoms = []
    for i, num in enumerate(numofatoms):
        # https://gitlab.com/ase/ase/blob/master/ase/io/vasp.py
        numofatoms[i] = int(num)
        [atoms.append(atomtypes[i]) for na in range(numofatoms[i])]

    if not vasp5:
        line_coortype = 6
    else:
        line_coortype = 7

    # TODO: Supporting Cartesian coordinates vasp input
    coortype = lines[line_coortype].split()[0]
    if coortype[0] in "sS":
        warnings.warn("Sorry! Selective dynamics "
                      "are not supported now", FutureWarning, stacklevel=2)
        line_coortype += 1
        coortype = lines[line_coortype].split()[0]

    if coortype[0] in "cCkK":
        line_first_pos = line_coortype + 1
        iscart=True
    else:
        iscart =False

    if coortype[0] in "dD":
        line_first_pos = line_coortype + 1

    positions = []
    total_atoms = sum(numofatoms)
    for i in range(line_first_pos, line_first_pos + total_atoms):
        s = lines[i].split()
        vec = float(s[0]), float(s[1]), float(s[2])
        positions.append(vec)
    if iscart:
        positions = numpy.dot(numpy.array(positions),numpy.linalg.inv(lattice))
    return Cell(lattice, positions, atoms)


def write_vasp(cell, filename='POSCAR', suffix='.vasp', long_format=True):
    """
    write vasp POSCAR type into file, vasp5 format only.
    always write atoms sorted POSCAR.

    parameters:

    cell: Cell object, the Cell that you wanna write into vasp POSCAR.
    filename: string, filename of output file, default='POSCAR'
    suffix: string, suffix of filename, default='.vasp'
    long_format: bool, if True format %.16f will be write, else %.6f
    ref: https://gitlab.com/ase/ase/blob/master/ase/io/vasp.py

    if optional parameters (filename and suffix) are not set,
    the filename will be 'POSCAR.vasp'
    """
    # TODO: write Cartesian coor POSCAR
    filname_suffix = ''.join([filename, suffix])
    with open(filname_suffix, "w") as f:
        f.write(_write_string(cell, long_format))


def _write_string(cell, long_format, print_vacc=False):
    """
    _write_string make io easy to be tested.

    return: string represtent POSCAR
    """

    # 对原子种类合并排序，用以产生体系名称和原子顺序数目和正确的坐标排序
    # sorted is a list of tuple(atom, na)
    atoms_dict = collections.Counter(cell.atoms)
    if not print_vacc:
        del atoms_dict[0]
    sorted_symbols = sorted(atoms_dict.items(), key=operator.itemgetter(0))

    list_symbols = ["{:}{:}".format(get_symbol(atom), na)
                    for atom, na in sorted_symbols]

    comment = ' '.join(list_symbols)
    comment += '\n'

    scale = '{:9.6f}'.format(1.0)
    scale += '\n'

    lattice_string = ""
    if long_format:
        latt_form = '21.16f'
    else:
        latt_form = '11.6f'

    for vec in cell.lattice:
        lattice_string += ' '
        for v in vec:
            lattice_string += '{:{form}}'.format(v, form=latt_form)
        lattice_string += '\n'

    # atom types and their numbers
    atom_types = ' '.join([get_symbol(i[0]) for i in sorted_symbols])
    atom_types += '\n'

    atom_numbers = ' '.join([str(i[1]) for i in sorted_symbols])
    atom_numbers += '\n'

    # TODO: write Cartesian coor
    coor_type = 'Direct\n'

    # argsort atoms and resort coor
    idx = numpy.argsort(cell.atoms)
    coord = cell.positions[idx]
    atoms = cell.atoms[idx]
    positions_string = ""
    if long_format:
        pos_form = '19.16f'
    else:
        pos_form = '9.6f'

    for i, vec in enumerate(coord):
        if atoms[i] == 0:
            continue
        positions_string += ' '
        for v in vec:
            positions_string += '{:{form}}'.format(v, form=pos_form)
        positions_string += ' ' + get_symbol(atoms[i])
        positions_string += '\n'

    poscar_string = ''.join([comment,
                             scale,
                             lattice_string,
                             atom_types,
                             atom_numbers,
                             coor_type,
                             positions_string])
    return poscar_string
