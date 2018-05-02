# coding: utf-8
# Distributed under the terms of the MIT License.
from .model import AppModel
from ababe.stru.element import Specie
from ababe.stru.scaffold import GeneralCell
from ababe.stru.sogen import OccupyGenerator
from ababe.io.io import GeneralIO
from ababe.stru.restriction import MinDistanceRestriction

import numpy as np
import os

class App(AppModel):

    def __init__(self, infile, comment, element, speckle, trs,
                 refined, outmode):
        # read comment & zoom from setting file first
        # if not exist, read from cmd args, then default
        self.cell = GeneralIO.from_file(infile)

        self.comment = comment or self.cell.comment

        self.element = element

        # Get number and index of target element
        num = self.cell.numbers
        if element is None:
            tgt_ele = int(num[1])
        else:
            tgt_ele = Specie(element).Z
        tgt_ele_index = np.where(num == tgt_ele)[0]
        self.n = tgt_ele_index.size

        self.s1, self.s2 = speckle
        # self.ele for function all-speckle-gen-of-ele in run
        self.ele = Specie.from_num(tgt_ele)

        # if there no restriction given then no restriction
        if trs != ():
            self.tr = trs[0]
        else:
            self.tr = None

        self.refined = refined
        self.outmode = outmode

    def run(self):
        # Create directory contain POSCARs
        import random
        import string

        rd_suffix = ''.join(random.choices(string.ascii_uppercase
                                           + string.digits, k=5))
        working_path = os.getcwd()
        out_dir = os.path.join(working_path,
                                   'STRUCTURES_{0:}_{1:}'.format(self.comment,
                                                              rd_suffix))
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        else:
            shutil.rmtree(out_dir)
            os.makedirs(out_dir)

        ogg = OccupyGenerator(self.cell)

        if self.tr is not None:
            tr = (Specie(self.tr[0]), self.tr[1])
            applied_restriction = MinDistanceRestriction(tr)


        for n1 in range(1, self.n - 1):
            for n2 in range(1, self.n - n1):
                g = ogg.gen_nodup_trinary_alloy(Specie(self.s1), n1,
                                                Specie(self.s2), n2)

                for n_count, c in enumerate(g):
                    if self.tr is not None:
                        condition = c.is_primitive() and applied_restriction.is_satisfied(c)
                    else:
                         condition = c.is_primitive()

                    if condition:
                        if self.refined:
                            c = c.get_refined_pcell()
                        out = GeneralIO(c)

                        f_suffix = ''.join(random.choices(string.ascii_uppercase
                                                          + string.digits, k=4))
                        ofname = "STRUCTURE_{:}_{:}.{:}".format(c.comment, f_suffix, self.outmode)
                        lastpath = os.path.join(out_dir, ofname)
                        out.write_file(lastpath)
