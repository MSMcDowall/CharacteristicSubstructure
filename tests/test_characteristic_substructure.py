import source.characteristic_substructure as cs
import unittest

class CharacteristicSubstructureTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_data_input(self):
        # What happens when file is empty?
        # What happens when file isn't smiles?
        # Looping over smiles correctly
        pass

    def test_find_graph_paths(self):
        # Incorrect SMILES format?
        # update molecules and path variables correctly
        pass

    def test_find_representative_paths(self):
        # Check only correct length paths are located
        # Check search for matching path
        # Check maths
        pass

    def test_find_representative_structures(self):
        pass

    def test_find_characteristic_substructure(self):
        pass

    def test_add_structure_to_characteristic(self):
        pass

    def test_create_structure(self):
        pass

    def test_create_multiple_structures(self):
        pass

    def test_create_nx_graph(self):
        # What happens if there is no nx_structure in place?
        pass

    def test_nx_isomorphism(self):
        pass

    def test_check_structure_duplicates(self):
        pass

    def test_swap_molecule_vertices(self):
        pass

    def test_swap_path_structure(self):
        pass