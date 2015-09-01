from ...source import algorithm
import unittest
import os


class AlgorithmFeatureTestCase(unittest.TestCase):
    def setUp(self):
        self.tester = algorithm.CSAlgorithm("")

    def test_find_all_representative_structures(self):
        pass

    def test_find_characteristic_substructure(self):
        pass

    def test_argument_input(self):
        pass

    def test_result_from_empty_file(self):
        print os.path.relpath("test/test_data/EmptyFile.txt", "source")
        self.tester.smiles_file = os.path.relpath("test", "source")
        pass