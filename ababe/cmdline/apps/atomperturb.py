# coding: utf-8
# Distributed under the terms of the MIT License.
from .model import AppModel
from ababe.stru.scaffold import ModifiedCell
from ababe.stru.element import Specie
from ababe.io.io import GeneralIO

import os
import numpy as np

class App(AppModel):

    def __init__(self, infile, radius):
        gcell = GeneralIO.from_file(infile)
        self.infile = infile
        self.basefname = os.path.basename(infile)
        self.mcell = ModifiedCell.from_gcell(gcell)

        self.radius = radius

    def run(self):
        import tempfile
        working_path = os.getcwd()

        self.mcell.perturb(self.radius)
        gcell = self.mcell.to_gcell()

        out = GeneralIO(gcell)
        ofname = "{:}_PURB.vasp".format(self.basefname.split('.')[0])

        print("PROCESSING: {:}".format(self.infile))
        out.write_file(ofname)
