import unittest

from pyabc.crystal.structure import Cell
from pyabc.crystal.derive import configurations_nonredundant


class TestDerive(unittest.TestCase):

    def test_conf_non_redun(self):
        fcc_latt = [0, 5, 5,
                    5, 0, 5,
                    5, 5, 0]
        fcc_pos = [(0, 0, 0)]
        fcc_atoms = [0]
        fcc_pcell = Cell(fcc_latt, fcc_pos, fcc_atoms)
        configurations_nonredundant(fcc_pcell, [[3, 5]], 8)
        # sc_latt = [5, 0, 0,
        #            0, 5, 0,
        #            0, 0, 5]
        # sc_pos = [(0, 0, 0)]
        # sc_atoms = [0]
        # sc_pcell = Cell(sc_latt, sc_pos, sc_atoms)
        # configurations_nonredundant(sc_pcell, [[3, 5]], 4)


if __name__ == "__main__":
    import nose2
    nose2.main()
