from pyabc.crystal.utils import non_dup_hnfs


def non_dup_lattice(pcell, volume=1, symprec=1e-5, comprec=1e-5):
    """
    non_dup_lattice return all non duplicated hnf extend cells.

    parameters:

    pcell: Cell object, The primitive cell to be extended
    volume: int, Extend to how large supercellself, default=1
    symprec: int, symmetry precision
                When finding duplicated hnfs the precesion, default=1e-5
    comprec: float, compare precision
                When finding the rotations symmetry of primitive cell, defalut=1e-5

    return:

    A list of Cell objects.
    """
    return [pcell.extend(hnf) for hnf in non_dup_hnfs(pcell, volumn, symprec, comprec)]
