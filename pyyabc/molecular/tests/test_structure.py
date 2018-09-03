import unittest
import numpy

from pyyabc.molecular.structure import Molecular


class TestMol(unittest.TestCase):

    def test_init(self):
        positions = [0, 0, 0]
        atoms = ['NaN_10']
        mole = Molecular(positions, atoms)
        self.assertEqual(mole.atoms.tolist(), [1010])

    def test_get_permutations(self):
        pos = numpy.array([0.5288,  0.1610,  0.9359,
                           0.0000,  0.0000,  0.0000,
                           0.2051,  0.8240, -0.6786,
                           0.3345, -0.9314, -0.4496,
                           -1.0685, -0.0537, 0.1921]).reshape((-1, 3))
        atoms = ['H', 'C', 'H', 'H', 'H']
        methane = Molecular(pos, atoms)

        symprec = 1e-2
        perm = methane.get_symmetry_permutation(symprec)
        num_sym = numpy.shape(perm)[0]
        # methane have 12 symmetry operation
        expected_num_sym = 12
        self.assertEqual(num_sym, expected_num_sym)

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
