# coding: utf-8
# Distributed under the terms of the MIT License.
import numpy as np
import sys
from collections import Counter, OrderedDict
from operator import itemgetter

from ababe.stru.element import Specie
from ababe.stru.scaffold import GeneralCell

class VaspInput(object):

    def __init__(self, filename):
        with open(filename, "r") as f:
            content = f.read()

        self.latt, self.pos, self.numbers = self.from_string(content)

    @staticmethod
    def from_string(content):
        # move empty line
        lines = [l for l in content.split('\n') if l.rstrip()]

        comment = lines[0]
        zoom = float(lines[1])
        lattice = np.around(np.array([[float(i) for i in line.split()]
                                                     for line in lines[2:5]]),
                                     decimals=6)
        if zoom < 0:
            # In vasp, a negative scale factor is treated as a volume. We need
            # to translate this to a proper lattice vector scaling.
            vol = abs(np.linalg.det(lattice))
            lattice *= (-zoom / vol) ** (1 / 3)
        else:
            lattice *= zoom

        nsymbols = [Specie(s).Z for s in lines[5].split()]
        natoms = [int(i) for i in lines[6].split()]
        numbers = []
        for i, s in enumerate(natoms):
            numbers += s * [nsymbols[i]]
        numbers = np.array(numbers)
        nsite = sum(natoms)

        lpos = 8    # first line of coordinates initial
        postype = lines[7].split()[0]

        # Selective dynamics
        if postype[0] in "Ss":
            lpos += 1

        positions = np.around(np.array([[float(i) for i in line.split()[0:3]]
                                                     for line in lines[lpos:lpos+nsite+1]]),
                              decimals=6)

        return lattice, positions, numbers

    def get_cell(self):
        return GeneralCell(self.latt, self.pos, self.numbers)


class VaspOutput(object):

    def __init__(self, gcell, comment=None, zoom=1):
        self.lattice = gcell.lattice
        self.positions = gcell.positions
        self.numbers = gcell.numbers

        self.atoms_name_list = list(map(lambda x: Specie.to_name(x),
                                        list(self.numbers)))

        self.zoom = zoom
        self.comment = comment

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

    def write_file(self, filename):
        """
        Writes POSCAR to a file.
        """
        with open(filename, "w") as f:
            f.write(self.get_string())

    def __repr__(self):
        return self.get_string()

    def __str__(self):
        """
        String representation of Poscar file
        """
        # print(self.get_string())
        return self.get_string()
