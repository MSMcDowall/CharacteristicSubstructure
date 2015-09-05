from ...characteristic_substructure.smiles_parser import Parser
import unittest


class ParserTestCase(unittest.TestCase):
    """
    Feature (end-to-end) tests for the SMILES parser.
    A SMILES string is given as input and the received molecule is tested to ensure it has been created correctly.
    """
    def setUp(self):
        self.parser = Parser()

    def tearDown(self):
        self.parser = Parser()

    def test_organic_creation(self):
        self.m = self.parser.parse_smiles('CN')
        first = self.m.vertex_from_position(0)
        second = self.m.vertex_from_position(1)
        self.assertEqual(first.element, 'C')
        self.assertEqual(second.element, 'N')
        self.assertTrue(self.m.adjacency_dictionary[first][second].single)

    def test_square_attributes(self):
        self.m = self.parser.parse_smiles('[12C@H2--][N@@H+3]')
        first = self.m.vertex_from_position(0)
        second = self.m.vertex_from_position(1)
        self.assertEqual(first.element, 'C')
        self.assertEqual(first.isotope, '12')
        self.assertEqual(first.hydrogen, '2')
        self.assertEqual(first.charge, '-2')
        self.assertEqual(second.element, 'N')
        self.assertEqual(second.hydrogen, '1')
        self.assertEqual(second.charge, '+3')
        self.assertTrue(self.m.adjacency_dictionary[first][second].single)

    def test_aromatic_bonds(self):
        self.m = self.parser.parse_smiles('ccN')
        first = self.m.vertex_from_position(0)
        second = self.m.vertex_from_position(1)
        third = self.m.vertex_from_position(2)
        self.assertTrue(first.aromatic)
        self.assertTrue(second.aromatic)
        first_bond = self.m.contains_edge(first, second)
        self.assertTrue(first_bond.aromatic)
        second_bond = self.m.contains_edge(second, third)
        self.assertTrue(second_bond.single)

    def test_ring_bond(self):
        self.m = self.parser.parse_smiles('C1CCC1')
        first = self.m.vertex_from_position(0)
        last = self.m.vertex_from_position(3)
        self.assertTrue(self.m.contains_edge(first, last))
        self.assertTrue(self.m.adjacency_dictionary[first][last].single)

    def test_branch(self):
        self.m = self.parser.parse_smiles('C(OC)C')
        first = self.m.vertex_from_position(0)
        root = self.m.vertex_from_position(1)
        branch_end = self.m.vertex_from_position(2)
        last = self.m.vertex_from_position(3)
        self.assertTrue(self.m.contains_edge(first, root))
        self.assertTrue(self.m.contains_edge(first, last))
        self.assertTrue(self.m.contains_edge(root, branch_end))
        self.assertFalse(self.m.contains_edge(branch_end, last))

    def test_bond(self):
        self.m = self.parser.parse_smiles('C-C=C#C$C:C')
        first = self.m.vertex_from_position(0)
        second = self.m.vertex_from_position(1)
        third = self.m.vertex_from_position(2)
        fourth = self.m.vertex_from_position(3)
        fifth = self.m.vertex_from_position(4)
        sixth = self.m.vertex_from_position(5)
        self.assertTrue(self.m.contains_edge(first, second).single)
        self.assertTrue(self.m.contains_edge(second, third).double)
        self.assertTrue(self.m.contains_edge(third, fourth).triple)
        self.assertTrue(self.m.contains_edge(fourth, fifth).quadruple)
        self.assertTrue(self.m.contains_edge(fifth, sixth).aromatic)

    def test_dot(self):
        self.m = self.parser.parse_smiles('C.N')
        self.c = self.m.vertex_from_position(0)
        self.n = self.m.vertex_from_position(1)
        self.assertFalse(self.m.contains_edge(self.c, self.n))

if __name__ == '__main__':
    unittest.main()