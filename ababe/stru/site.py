# coding: utf-8
# Distributed under the terms of the MIT License.
from ababe.stru.element import Specie
import numpy as np

class Site(object):

    def __init__(self, position, ele):
        self._position = tuple(position)
        if isinstance(ele, Specie):
            self._element = ele
        elif isinstance(ele, int):
            self._element = Specie.from_num(ele)
        else:
            self._element = Specie(ele)

    def __eq__(self, other):
        return self.position == other.position and self.element == other.element

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, atuple):
        self._position = tuple(atuple)

    @property
    def element(self):
        return self._element

    @element.setter
    def element(self, ele):
        self._element = ele
