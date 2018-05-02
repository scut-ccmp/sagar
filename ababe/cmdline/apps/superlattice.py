# coding: utf-8
# Distributed under the terms of the MIT License.
from .model import AppModel
from ababe.io.io import GeneralIO
from ababe.stru.grid import SuperLatticeGenerator, SuperLatticeGenerator2D
import numpy as np
import os
import shutil

class App(AppModel):

    def __init__(self, infile, comment, volumn, ld, outmode):
        self.ucell = GeneralIO.from_file(infile)
        # check whether input a unit cell
        if not self.ucell.is_primitive():
            raise ValueError('Lattice in setting file are not primitive.\n'
                             'You can reinput OR get the primitive lattice\n'
                             'by using runati pcell <file>\n')
        self.v = volumn

        self.comment = comment or "default"

        self.outmode = outmode

        if ld:
            LatticeGen = SuperLatticeGenerator2D
        else:
            LatticeGen = SuperLatticeGenerator

        self.hnfs = LatticeGen.hnfs_from_n(self.ucell, self.v)

    def run(self):
        import random
        import string

        # Create a dir contains suplat files
        working_path = os.getcwd()
        suplat_dir = os.path.join(working_path,
                                  'SUPLAT_{:}'.format(self.comment))
        if not os.path.exists(suplat_dir):
            os.makedirs(suplat_dir)
        else:
            shutil.rmtree(suplat_dir)
            os.makedirs(suplat_dir)

        for hnf in self.hnfs:
            rd_suffix = ''.join(random.choices(string.ascii_uppercase
                                               + string.digits, k=4))
            sl_origin = hnf.to_general_cell()
            sl = sl_origin.get_shaped_cell()

            out = GeneralIO(sl)

            ofname = "SUPLAT_{:}_{:}.{:}".format(self.v, rd_suffix, self.outmode)
            lastpath = os.path.join(suplat_dir, ofname)
            out.write_file(lastpath)
