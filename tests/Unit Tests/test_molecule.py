import source.molecule as molecule
import unittest

class AtomTestCase(unittest.TestCase):
    def test_atom_creation(self):
        self.a = molecule.Atom('C', hydrogen=2)
        self.assertEqual(self.a.element, 'C')       # Test for access to base class methods
        self.assertEqual(self.a.hydrogen, 2)

    def test_aromatic_atom_creation(self):
        self.a = molecule.Atom('c', 12, aromatic=True)
        self.assertEqual(self.a.element, 'c')
        self.assertEqual(self.a.isotope, 12)

class BondTestCase(unittest.TestCase):
    def setUp(self):
        self.x = molecule.Atom('C')
        self.y = molecule.Atom('N')

    def test_single_bond_creation(self):
        self.b = molecule.Bond(self.x, self.y, single=True)
        self.assertEqual(self.b.endpoints, (self.x, self.y))       # Test for access to base class methods
        self.assertEqual(self.b.opposite(self.x), self.y)
        self.assertEqual(self.b.single, True)

    def test_double_bond_creation(self):
        self.b = molecule.Bond(self.x, self.y, double=True)
        self.assertEqual(self.b.endpoints, (self.x, self.y))
        self.assertEqual(self.b.opposite(self.x), self.y)

    def test_triple_bond_creation(self):
        self.b = molecule.Bond(self.x, self.y, triple=True)
        self.assertEqual(self.b.endpoints, (self.x, self.y))
        self.assertEqual(self.b.opposite(self.x), self.y)

    def test_quadruple_bond_creation(self):
        self.b = molecule.Bond(self.x, self.y, quadruple=True)
        self.assertEqual(self.b.endpoints, (self.x, self.y))
        self.assertEqual(self.b.opposite(self.x), self.y)

    def test_aromatic_bond_creation(self):
        self.b = molecule.Bond(self.x, self.y, aromatic=True)
        self.assertEqual(self.b.endpoints, (self.x, self.y))
        self.assertEqual(self.b.opposite(self.x), self.y)


class MoleculeTestCase(unittest.TestCase):
    def setUp(self):
        self.mole = molecule.Molecule('string')
        self.x = self.mole.add_atom('C')
        self.y = self.mole.add_atom('C')

    def test_molecule_creation(self):
        self.assertEqual(self.mole.size, 2)     # Test for access to base class methods

    def test_add_atom(self):
        self.assertEqual(self.mole.vertices()[0].element, 'C')

    def test_add_aromatic_atom(self):
        self.mole.add_aromatic_atom('c')
        self.assertEqual(self.mole.vertices()[2].element, 'c')

    def test_add_single_bond(self):
        self.e = self.mole.add_single_bond(self.x, self.y)
        self.assertEqual(self.mole.adjacency_dictionary[self.x][self.y], self.e)

    def test_add_double_bond(self):
        self.e = self.mole.add_double_bond(self.x, self.y)
        self.assertEqual(self.mole.adjacency_dictionary[self.x][self.y], self.e)

    def test_add_triple_bond(self):
        self.e = self.mole.add_triple_bond(self.x, self.y)
        self.assertEqual(self.mole.adjacency_dictionary[self.x][self.y], self.e)

    def test_add_quadruple_bond(self):
        self.e = self.mole.add_quadruple_bond(self.x, self.y)
        self.assertEqual(self.mole.adjacency_dictionary[self.x][self.y], self.e)

    def test_add_aromatic_bond(self):
        self.aromatic_first = self.mole.add_aromatic_atom('c')
        self.aromatic_second = self.mole.add_aromatic_atom('n')
        self.e = self.mole.add_aromatic_bond(self.aromatic_first, self.aromatic_second)
        self.assertEqual(self.mole.adjacency_dictionary[self.aromatic_first][self.aromatic_second], self.e)

if __name__ == '__main__':
    unittest.main()