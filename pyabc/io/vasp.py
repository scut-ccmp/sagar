import numpy
import warnings

from pyabc.crystal.structure import Cell


def read(filename='POSCAR'):
    """
    Import POSCAR/CONTCAR or filename with .vasp suffix

    parameter:
    filename: string, the filename

    return: Cell object.
    """
    # TODO: read velocities, now not supported.
    with open(filename, "r") as f:
        # _read_string return Cell object.
        return _read_cell_from_string(f.read())


def _read_cell_from_string(data):
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
        warnings.warn("symbols of elements in fifth line are missing,"
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
        warnings.warn("Sorry! Selective dynamics"
                      "are not supported now", FutureWarning, stacklevel=2)
        line_coortype += 1
    if coortype[0] in "cCkK":
        raise ValueError("Sorry! Cartesian coordinates"
                         "are not supported now,"
                         "modify your input file.")
    if coortype[0] in "dD":
        line_first_pos = line_coortype + 1

    positions = []
    total_atoms = sum(numofatoms)
    for i in range(line_first_pos, line_first_pos + total_atoms):
        s = lines[i].split()
        vec = float(s[0]), float(s[1]), float(s[2])
        positions.append(vec)

    return Cell(lattice, positions, atoms)


def write():
    pass
