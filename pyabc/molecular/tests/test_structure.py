import unittest
import numpy
from pyabc.molecular.structure import Molecular


class TestMol(unittest.TestCase):

    # def setUp(self):
    #     # positions = numpy.loadtxt('z12.txt')
    #     # atoms = ['C'] * 20
    #     # self.molecular = Molecular(positions, atoms)

    def test_init(self):

        positions = [0, 0, 0]
        atoms = ['NaN_10']
        mole = Molecular(positions, atoms)
        self.assertEqual(mole.atoms.tolist(), [1010])

    def test_get_permutations(self):
        pos = numpy.loadtxt('pos.txt')
        atoms = ['C'] * 60
        c60 = Molecular(pos, atoms)
        self.assertEqual(numpy.shape(c60.get_symmetry_permutation(pres=0.05))
                         [0], 120)

    def test_check(self):
        # c-c键和c-o键均距离不大于0.1A
        cco_pos = [2.0, 0.0, 0.0,
                   -2.0, 0.0, 0.0,
                   0.0, 0.0, 0.0]
        cco_atoms = [6, 6, 8]
        mol = Molecular(cco_pos, cco_atoms)
        self.assertTrue(mol.check(elements=None, limit=0.1))

        # 测试元素c-c距离过近, 但只测o元素时候不会过近
        cco_pos = [2.0, 0.0, 0.0,
                   2.01, 0.0, 0.0,
                   0.0, 0.0, 0.0]
        cco_atoms = [6, 6, 8]
        mol = Molecular(cco_pos, cco_atoms)
        self.assertFalse(mol.check(elements=['C'], limit=0.1))
        self.assertTrue(mol.check(elements=['O'], limit=0.1))


if __name__ == "__main__":
    import nose2
    nose2.main()
