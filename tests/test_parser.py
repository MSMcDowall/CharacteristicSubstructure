import source.parser as parser
import source.molecule as molecule
import unittest

class TestParser(unittest.TestCase):
    def setUp(self):
        self.m = molecule.Molecule('')

    def tearDown(self):
        self.m.clear()

    def test_organic(self):
        self.m = parser.parse_smiles('CN')
        first = self.m.vertices[0]
        second = self.m.vertices[1]
        self.assertEqual(first.element, 'C')
        self.assertEqual(second.element, 'N')
        self.assertIsInstance(self.m._vertices[first][second], molecule.SingleBond)

    # Through debugger does not enter into add_single_bond but when run normally it does
    """
    def test_square(self):
        #square = '[N@@H+3]'
        self.m = parser.parse_smiles('[12C@H2--]')
        first = self.m.vertices[0]
        #second = self.m.vertices[1]
        self.assertEqual(first.element, 'C')
        self.assertEqual(first.isotope, '12')
        self.assertEqual(first.hydrogen, '2')
        self.assertEqual(first.charge, '-2')
    """

    def test_ring(self):
        pass

    def test_branch_start(self):
        pass

    def test_branch_end(self):
        pass

    def test_bond(self):
        pass

    def test_dot(self):
        self.m = parser.parse_smiles('C.N')
        pass


def main():
    unittest.main()


if __name__ == '__main__':
    main()