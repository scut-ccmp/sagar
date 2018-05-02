# coding: utf-8
# Distributed under the terms of the MIT License.

from ababe.stru.element import Specie
from collections import Counter, OrderedDict
from operator import itemgetter

import ruamel.yaml as yaml
# import yaml


class YamlOutput(object):

    def __init__(self, gcell, zoom=1):
        self.lattice = gcell.lattice
        self.positions = gcell.positions
        self.numbers = gcell.numbers

        self.atoms_name_list = list(map(lambda x: Specie.to_name(x),
                                        list(self.numbers)))

        self.zoom = zoom

    def get_string(self):
        """
        Returns:
            String representation of YAML type
        """
        d = Counter(self.atoms_name_list)
        ordered_atoms = OrderedDict(sorted(d.items(),
                                           key=lambda x: Specie(x[0]).Z))
        if 'G' in ordered_atoms:
            del ordered_atoms['G']

        comment = ''.join(['{}{}'.format(k, v)
                           for k, v in ordered_atoms.items()])

        output = {'comment':    comment,
                  'lattice':    self.lattice.tolist(),
                  'positions':  self.positions.tolist(),
                  'numbers':    self.numbers.tolist(),
                  'zoom':       self.zoom}

        return yaml.dump(output)

    def write(self, filename):
        """
        Writes YAML to a file.
        """
        with open(filename, "w") as f:
            f.write(self.get_string())

    def __repr__(self):
        return self.get_string()

    def __str__(self):
        """
        String representation of structure yaml file
        """
        # print(self.get_string())
        return self.get_string()


class VaspPOSCAR(object):

    def __init__(self, gcell, zoom=1):
        self.lattice = gcell.lattice
        self.positions = gcell.positions
        self.numbers = gcell.numbers

        self.atoms_name_list = list(map(lambda x: Specie.to_name(x),
                                        list(self.numbers)))

        self.zoom = zoom

    def get_string(self, direct=True):
        """
        Returns:
            String representation of POSCAR.
        """

        # latt = self.lattice
        d = Counter(self.atoms_name_list)
        ordered_atoms = OrderedDict(sorted(d.items(),
                                           key=lambda x: Specie(x[0]).Z))
        if 'G' in ordered_atoms:
            del ordered_atoms['G']

        comment = ''.join(['{}{}'.format(k, v)
                           for k, v in ordered_atoms.items()])

        lines = [comment, str(self.zoom)]
        # lattice string
        for c in self.lattice:
            line = " ".join("{:10.6f}".format(p) for p in c)
            lines.append(line)

        lines.append(" ".join([str(x) for x in ordered_atoms.keys()]))
        lines.append(" ".join([str(x) for x in ordered_atoms.values()]))

        zipped_list = list(zip(self.numbers, self.positions,
                               self.atoms_name_list))
        sorted_position = sorted(zipped_list, key=itemgetter(0))

        lines.append("direct" if direct else "cartesian")

        # sort the positions by atoms seqence
        for (i, pos, site) in sorted_position:
            if not i == 0:
                line = " ".join(["{:10.6f}".format(p) for p in pos])
                line += " " + site
                lines.append(line)

        return "\n".join(lines) + "\n"

    def __repr__(self):
        return self.get_string()

    def __str__(self):
        """
        String representation of Poscar file
        """
        # print(self.get_string())
        return self.get_string()

    def write(self, filename):
        """
        Writes POSCAR to a file.
        """
        with open(filename, "w") as f:
            f.write(self.get_string())


class GeneralCellInput(object):
    """ read files and get GeneralCell

        Input: filename
    """

    def __init__(self, filename):
        pass

    def get_gcell(self):
        pass
