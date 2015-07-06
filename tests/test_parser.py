import source.parser as parser
import source.molecule as molecule
import unittest

class ParserTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_organic(self):
        self.m = parser.parse_smiles('CN')
        first = self.m.vertices[0]
        second = self.m.vertices[1]
        self.assertEqual(first.element, 'C')
        self.assertEqual(second.element, 'N')
        self.assertIsInstance(self.m._vertices[first][second], molecule.SingleBond)

    def test_square(self):
        self.m = parser.parse_smiles('[12C@H2--][N@@H+3]')
        first = self.m.vertices[0]
        second = self.m.vertices[1]
        self.assertEqual(first.element, 'C')
        self.assertEqual(first.isotope, '12')
        self.assertEqual(first.hydrogen, '2')
        self.assertEqual(first.charge, '-2')
        self.assertEqual(second.element, 'N')
        self.assertEqual(second.hydrogen, '1')
        self.assertEqual(second.charge, '+3')
        self.assertIsInstance(self.m._vertices[first][second], molecule.SingleBond)

    def test_ring(self):
        self.m = parser.parse_smiles('c1ccccc1')
        first = self.m.vertices[0]
        last = self.m.vertices[5]
        self.assertIsInstance(self.m._vertices[first][last], molecule.AromaticBond)

    def test_branch(self):
        self.m = parser.parse_smiles('C(OC)C')



    def test_bond(self):
        pass

    def test_dot(self):
        self.m = parser.parse_smiles('C.N')
        self.c = self.m.vertices[0]
        self.n = self.m.vertices[1]
        self.assertFalse(self.m.contains_edge(self.c, self.n))


if __name__ == '__main__':
    unittest.main()