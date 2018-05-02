import numpy

class Cell(object):
    """
    Cell object represent a crystal structure.

    parameters:

    lattice: 3x3 1D-list, lattice of cell

    positions: n tuples(x,y,z) in fraction.

    atoms: list of numbers, represent atom in periodic table.
    """
    def __init__(self, lattice, positions, atoms):
        self._lattice = numpy.array(lattice).reshape((3,3))
        self._numbers = len(positions)
        self._positions = numpy.array(positions)
        self._atoms = numpy.array(atoms)

    def extend(self, mat):
        latt = numpy.matmul()
        return self.__class__(lattice, positions, atoms)
