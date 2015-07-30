from source.smiles_parser import Parser
import source.molecule as molecule
import unittest

class ParserTestCase(unittest.TestCase):
    # TODO Currently focused on functionality testing; add unit tests
    def setUp(self):
        self.parser = Parser()

    def tearDown(self):
        self.parser = Parser()

    def test_organic(self):
        self.m = self.parser.parse_smiles('CN')
        first = self.m.vertices[0]
        second = self.m.vertices[1]
        self.assertEqual(first.element, 'C')
        self.assertEqual(second.element, 'N')
        self.assertIsInstance(self.m._vertices[first][second], molecule.SingleBond)

    def test_square(self):
        self.m = self.parser.parse_smiles('[12C@H2--][N@@H+3]')
        first = self.m.vertices[0]
        second = self.m.vertices[1]
        self.assertEqual(first.element, 'C')
        self.assertEqual(first.isotope, '12')
        self.assertEqual(first.hydrogen, '2')
        self.assertEqual(first.charge, '-2')
        self.assertEqual(second.element, 'N')
        self.assertEqual(second.hydrogen, '1')
        self.assertEqual(second.charge, '+3')
        print self.m._vertices[first]
        self.assertIsInstance(self.m._vertices[first][second], molecule.SingleBond)

    def test_aromatic(self):
        self.m = self.parser.parse_smiles('ccN')
        # Aromatic atoms 1 and 2
        # Aromatic bond between 1 and 2
        # Single bond between 2 and 3

    # FIXME Faulty test, occasionally passes and fails without code change
    """
    def test_ring_bond(self):
        self.m = self.parser.parse_smiles('C1CCC1')
        first = self.m.vertices[0]
        last = self.m.vertices[3]
        print repr(first)
        print repr(last)
        print self.m._vertices[first]
        self.assertIsInstance(self.m._vertices[first][last], molecule.SingleBond)
    """

    def test_branch(self):
        self.m = self.parser.parse_smiles('C(OC)C')
        pass

    def test_bond(self):
        pass

    def test_dot(self):
        self.m = self.parser.parse_smiles('C.N')
        self.c = self.m.vertices[0]
        self.n = self.m.vertices[1]
        self.assertFalse(self.m.contains_edge(self.c, self.n))


if __name__ == '__main__':
    unittest.main()