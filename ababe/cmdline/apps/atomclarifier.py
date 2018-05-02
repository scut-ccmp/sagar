# coding: utf-8
# Distributed under the terms of the MIT License.
from .model import AppModel
from ababe.stru.scaffold import ModifiedCell
from ababe.stru.clarifier import VerboseAtomRemoveClarifier
from ababe.stru.element import Specie
from ababe.io.io import GeneralIO

import os
import numpy as np

class App(AppModel):

    def __init__(self, files, cenele, radius, ele, refined):
        self.files = files

        self.clarifier = VerboseAtomRemoveClarifier(Specie(cenele), radius, Specie(ele))
        self.refined = refined

    def run(self):
        import tempfile
        working_path = os.getcwd()

        for infile in self.files:
            basefname = os.path.basename(infile)

            # read
            gcell = GeneralIO.from_file(infile)

            mcell = ModifiedCell.from_gcell(gcell)
            new_mcell = self.clarifier.clarify(mcell)
            gcell = new_mcell.to_gcell()
            # todo: add feature- to convcell.
            if self.refined:
                gcell = gcell.get_refined_pcell()

            out = GeneralIO(gcell)
            ofname = "{:}_ACLR.vasp".format(basefname.split('.')[0])

            print("PROCESSING: {:}".format(infile))
            out.write_file(ofname)
