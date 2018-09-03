import numpy
from itertools import combinations

from pyabc.molecular.structure import Molecular
from pyabc.utils.core import _remove_redundant

class ConfigurationGenerator(object):
    '''
    这个类用于产生各种替换原子的需求
    '''

    def __init__(self, mol, symprec=1e-3):
        if not isinstance(mol, Molecular):
            raise TypeError(
                "want pyabc.molecular.structure.Molecular, got {:}".
                format(type(mol)))
        self.mol = mol
        self.perms = mol.get_symmetry_permutation(symprec)

    def get_configurations(self, sites, e_num):
        '''
        get_configurations output specific molecular
        for specific concentration.

        parameters:
        sites: list of (lists or tuples), represent element disorder of each sites
        e_num: tuple, number of atoms in disorderd sites.

        yield:

        a tuple
        tuple[0]: Cell object, a list of non-redundant configurations of certain volume supercell.
        tuple[1]: int object, degeneracy of the configuration in all configurations of this volume.
        '''
        perms = self.perms
        mol = self.mol
        for pa, d in _remove_redundant(mol, sites, perms, e_num=e_num):
            m = Molecular(pa[0], pa[1])
            yield (m, d)
