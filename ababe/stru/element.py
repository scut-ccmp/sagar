# coding: utf-8
# Distributed under the terms of the MIT License.
import os
import json
import pdb

from aenum import Enum

# Loads periodic_table database from json file pt.json
with open(os.path.join(os.path.dirname(__file__),
                            "pt.json"), "rt") as f:
    _pt_db = json.load(f)

class Specie(object):

    def __init__(self, symbol):
        self.symbol = "%s" % symbol
        d = _pt_db[symbol]

        self._Z = d["Atomic no"]
        self.atom_mass = d["Atomic mass"]
        self.atom_radius = d["Atomic radius"]

    @property
    def Z(self):
        return self._Z

    @property
    def name(self):
        return self.symbol

    @staticmethod
    def from_num(an):
        if an == 0:
            sp = GhostSpecie()
        else:
            key = [k for k, v in _pt_db.items() if v["Atomic no"] == an][0]
            sp = Specie(key)
        return sp

    @staticmethod
    def to_name(an):
        if an == 0:
            name = "G"
        else:
            key = [k for k, v in _pt_db.items() if v["Atomic no"] == an][0]
            name = key

        return name

    def __eq__(self, other):
        return isinstance(other, Specie) and self.Z == other.Z

    def __ne__(self, other):
        return not self.__eq__(other)


class GhostSpecie(Specie):
    """
    A special species for representing vacancies
    """

    def __init__(self):
        self.symbol = "G"
        self._Z = 0
        self.atom_mass = 0
        self.atom_radius = 0
