from ...source import algorithm
import unittest


class AlgorithmFeatureTestCase(unittest.TestCase):
    def setUp(self):
        self.tester = algorithm.CSAlgorithm()

    def test_result_from_empty_set(self):
        """
        If no SMILES are provided then no characteristic substructure is made
        :return:
        """
        smiles_set = [""]
        paths = self.tester.find_graphs_paths(smiles_set)
        cs = self.tester.find_characteristic_substructure(paths)
        self.assertTrue(cs.adjacency_dictionary == {})

    def test_result_from_single_common_chain(self):
        """
        Chain of carbons is common to all but the last structure which contains a carbon ring
        :return: None
        """
        smiles_set = ["CCCCCCNO", "CCCCCCPS", "CCCCCCNS", "CCCCCCO", "C1CCCCC1O"]
        paths = self.tester.find_graphs_paths(smiles_set)
        cs = self.tester.find_characteristic_substructure(paths)
        print cs.adjacency_dictionary_display()
        self.assertTrue(len(cs.adjacency_dictionary.keys()) == 6)
        middle = 0
        end = 0
        for vertex in cs.adjacency_dictionary:
            if len(cs.adjacency_dictionary[vertex]) == 2:
                middle += 1
            elif len(cs.adjacency_dictionary[vertex]) == 1:
                end += 1
        self.assertTrue(middle == 4)
        self.assertTrue(end == 2)

    def test_result_from_single_common_ring(self):
        """
        Carbon ring is common to most structures, ensure the structure is made correctly
        :return: None
        """
        smiles_set = ["C1CCCCC1NO", "C1CCCCC1PS", "C1CCCCC1NS", "C1CCCCC1O", "CCCCCC"]
        paths = self.tester.find_graphs_paths(smiles_set)
        cs = self.tester.find_characteristic_substructure(paths)
        print cs.adjacency_dictionary_display()
        self.assertTrue(len(cs.adjacency_dictionary.keys()) == 6)
        middle = 0
        end = 0
        for vertex in cs.adjacency_dictionary:
            if len(cs.adjacency_dictionary[vertex]) == 2:
                middle += 1
            elif len(cs.adjacency_dictionary[vertex]) == 1:
                end += 1
        self.assertTrue(middle == 6)
        self.assertTrue(end == 0)

    def test_result_from_repeated_common_structure_two_occurrences_disjoint(self):
        """
        SNO appears multiple times in most of the structures, either twice or three times. Always disconnected
        Ensure that no bond is made between the repeated structures, and that the CS contains two copies of SNO
        :return: None
        """
        self.tester.length_end = 3
        smiles_set = ["SNOCCCSNO", "SNOPPSNO", "SNOClSNOClClSNO", "SNOIIISNO", "CNO"]
        paths = self.tester.find_graphs_paths(smiles_set)
        cs = self.tester.find_characteristic_substructure(paths)
        print cs.adjacency_dictionary
        s = 0
        n = 0
        o = 0
        for vertex in cs.adjacency_dictionary:
            if vertex.label == 'S':
                s += 1
                neighbour = cs.adjacency_dictionary[vertex].keys()[0]
                self.assertTrue(neighbour.label == 'N')
                self.assertTrue(len(cs.adjacency_dictionary[vertex]) == 1)
            elif vertex.label == 'N':
                n += 1
                self.assertTrue(len(cs.adjacency_dictionary[vertex]) == 2)
            elif vertex.label == 'O':
                o += 1
                neighbour = cs.adjacency_dictionary[vertex].keys()[0]
                self.assertTrue(neighbour.label == 'N')
                self.assertTrue(len(cs.adjacency_dictionary[vertex]) == 1)
        self.assertTrue(s == 2)
        self.assertTrue(n == 2)
        self.assertTrue(o == 2)