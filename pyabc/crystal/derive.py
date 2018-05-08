from pyabc.crystal.utils import non_dup_hnfs


def cells_nonredundant(pcell, volume=1, symprec=1e-5, comprec=1e-5):
    """
    cells_nonredundant return all non duplicated hnf extend cells.

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

def configurations_nonredundant(pcell, volume=1, symprec=1e-5, comprec=1e-5):
    """
    在较大体积的构型中，可能有结构已经在小的体积中出现过。
    """
    pass
