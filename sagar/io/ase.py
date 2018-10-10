# -*- coding: utf-8 -*-
from ase import Atoms

from sagar.crystal.structure import Cell


def read(ase_atoms):
    if isinstance(ase_atoms, Atoms):
        lattice = ase_atoms.get_cell(complete=False)
        positions = ase_atoms.get_positions(wrap=False)
        atoms = ase_atoms.get_atomic_numbers()
        return Cell(lattice, positions, atoms)
    else:
        raise TypeError("The value is not an ase.Atoms object")


def write(cell):
    if isinstance(cell, Cell):
        ase_cell = cell.lattice
        ase_positions = cell.positions
        ase_numbers = cell.atoms
        return Atoms(cell=ase_cell,
                     scaled_positions=ase_positions,
                     numbers=ase_numbers,
                     pbc=(True, True, True))
    else:
        raise TypeError(
            "The value is not an sagar.crystal.structure Cell object")
