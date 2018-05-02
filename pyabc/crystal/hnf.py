import numpy

from math import sqrt
from itertools import product


def _factor(n):
    """
    https://rosettacode.org/wiki/Factors_of_an_integer#Python

    realize that factors come in pairs, the smaller of which is
    no bigger than sqrt(n)
    """
    factors = set()
    for x in range(1, int(sqrt(n)) + 1):
        if n % x == 0:
            factors.add(x)
            factors.add(n // x)
    return sorted(factors)

def _hnfs(det):
    for a in _factor(det):
        for d in _factor(det // a):
            f = det // a // d

            for b, c, e in product(range(d), range(f), range(f)):
                yield numpy.array([[a, b, c],
                                   [0, d, e],
                                   [0, 0, f]])

def hnf_cells(pcell, volume=1, symprec=1e-5):
    nodup_cells = []
    for hnf in _hnfs(volume):
        npc = pcell.extend(hnf) # or? extend(pcell, hnf)
        if npc not in nodup_cells:
            nodup_cells.append(npc)

    return nodup_cells
