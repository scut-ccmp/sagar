# coding: utf-8
# Distributed under the terms of the MIT License.

from ababe.stru.element import Specie
from collections import Counter, OrderedDict
from operator import itemgetter
from fnmatch import fnmatch


class GeneralIO(object):

    def __init__(self, gcell):
        self.gcell = gcell

    def write_file(self, fname=None, fmt=None):
        from ababe.io.vasp import VaspOutput
        from ababe.io.yaml import YamlOutput
        fname = fname or ""
        fmt = "" if fmt is None else fmt.lower()

        if fmt == "yaml" or fnmatch(fname.lower(), "*.yaml*"):
            writer = YamlOutput(self.gcell)
        elif fmt == "vasp" or fnmatch(fname.lower(), "*.vasp") or fnmatch(fname, "*POSCAR*"):
            writer = VaspOutput(self.gcell)

        if fname:
            writer.write_file(fname)
        else:
            print(writer.__str__())

    @classmethod
    def from_file(cls, fname, fmt=None):
        from ababe.io.vasp import VaspInput
        from ababe.io.yaml import YamlInput
        fname = fname or ""
        fmt = "" if fmt is None else fmt.lower()

        if fmt == "yaml" or fnmatch(fname.lower(), "*.yaml*"):
            reader = YamlInput(fname)
        elif fmt == "vasp" or fnmatch(fname.lower(), "*.vasp") or fnmatch(fname, "*POSCAR*"):
            reader = VaspInput(fname)

        return reader.get_cell()
