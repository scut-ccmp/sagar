# coding: utf-8
# Distributed under the terms of the MIT License.
import numpy as np

from ababe.stru.scaffold import ModifiedCell

class Clarifier(object):

    def __init__(self):
        pass

    def clarifier(self, modecell):
        pass


class AtomRemoveClarifier(Clarifier):

    def __init__(self, centers, r, element=None):
        self.centers = np.array(centers)
        self.element = element
        self.r = r

    def clarify(self, modcell):
        if not isinstance(modcell, ModifiedCell):
            raise ValueError('Please input a ModifiedCell')

        for center in self.centers:
            dsites = modcell.get_points_incell_insphere(center, self.r, self.element)
            modcell.remove_sites(list(dsites.keys()))
        return modcell


class VerboseAtomRemoveClarifier(Clarifier):

    def __init__(self, etobesub, r, element=None):
        self.etobesub = etobesub
        self.element = element
        self.r = r

    def clarify(self, modcell):
        if not isinstance(modcell, ModifiedCell):
            raise ValueError('Please input a ModifiedCell')

        for site in modcell:
            if site.element == self.etobesub:
                center = site.position
                dsites = modcell.get_points_incell_insphere(center, self.r, self.element)
                modcell.remove_sites(list(dsites.keys()))
        return modcell
