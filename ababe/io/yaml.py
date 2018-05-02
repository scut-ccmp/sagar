# coding: utf-8
# Distributed under the terms of the MIT License.
import numpy as np
import sys
from ruamel.yaml import YAML
from ruamel.yaml.compat import StringIO
from collections import Counter, OrderedDict

from ababe.stru.element import Specie
from ababe.stru.scaffold import GeneralCell

class MyYAML(YAML):
    def dump(self, data, stream=None, **kw):
        inefficient = False
        if stream is None:
            inefficient = True
            stream = StringIO()
        YAML.dump(self, data, stream, **kw)
        if inefficient:
            return stream.getvalue()

class YamlInput(object):

    def __init__(self, filename):
        with open(filename, "r") as f:
            content = f.read()

        self.latt, self.pos, self.numbers = self.from_string(content)

    @staticmethod
    def from_string(content):
        """
        The content is a string read from file and convert to
        this function. To get the lattice positions and numbers.
        """
        yaml = YAML(typ='safe')
        dicty = yaml.load(content)
        latt = np.around(np.array(dicty['lattice']) * dicty['zoom'], decimals=6)
        pos = np.around(np.array(dicty['positions']), decimals=6)
        numbers = np.array(dicty['numbers'])

        return latt, pos, numbers

    def get_cell(self):
        return GeneralCell(self.latt, self.pos, self.numbers)

class YamlOutput(object):

    def __init__(self, gcell):
        self.lattice = np.around(gcell.lattice, decimals=6)
        self.positions = np.around(gcell.positions, decimals=6)
        self.numbers = gcell.numbers

        atoms_name_list = list(map(lambda x: Specie.to_name(x),
                                   list(self.numbers)))
        d = Counter(atoms_name_list)
        ordered_atoms = OrderedDict(sorted(d.items(),
                                           key=lambda x: Specie(x[0]).Z))
        # remove Ghostatoms
        if 'G' in ordered_atoms:
            del ordered_atoms['G']

        self.comment = ''.join(['{}{}'.format(k, v)
                               for k, v in ordered_atoms.items()])

    def get_string(self):
        """
        Returns:
            String representation of yaml type
        """
        output = {'comment':    self.comment,
                  'lattice':    self.lattice.tolist(),
                  'positions':  self.positions.tolist(),
                  'numbers':    self.numbers.tolist(),
                  'zoom':       1}

        yaml = MyYAML(typ='safe')
        return yaml.dump(output)

    def write_file(self, filename):
        """
        write yaml to file
        """
        with open(filename, "w") as f:
            f.write(self.get_string())

    def __repr__(self):
        return self.get_string()

    def __str__(self):
        """
        String representation of structure yaml file
        """
        return self.get_string()
