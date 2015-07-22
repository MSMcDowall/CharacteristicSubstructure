# coding=utf-8
from parser import Parser

# Implementation of Finding Characteristic Substructures for Metabolite Classes
# Ludwig, Hufsky, Elshamy, Böcker

class CharacteristicSubstructure(object):
    def __init__(self, length_start=20, length_end=5, step=4, path_threshold=0.8, structure_threshold=0.8):
        self.length_start = length_start
        self.length_end = length_end
        self.step = step
        self.path_threshold = path_threshold
        self.structure_threshold = structure_threshold
        self.molecules = []                 # All the given molecules
        self.paths = {}                     # All the paths from all the molecules

    def find_graphs_paths(self, smiles_set):
        for smiles in smiles_set:
            molecule = Parser().parse_smiles(smiles)
            self.molecules.append(molecule)
            path_dict = molecule.find_all_paths()
            self.paths.update(path_dict)
        return self.paths

    def find_representative_paths(self, length):
        representative_paths = []
        for path in self.paths:
            if self.paths[path] == length:
                counter = 0
                for mole in self.molecules:
                    if path in mole.paths:
                        counter += 1
                if float(counter)/len(self.molecules) >= self.path_threshold:
                    representative_paths.append(path)
        return representative_paths

    def find_representative_structures(self, rep_paths):
        # Create structures from path
        # Test for subgraph in each molecule
        # Store location of each substructure that is isomorphic to representative structure as:
        # Dictionary within dictionary {path_structure: {molecule : position}}
        # Consider here what happens if there is more than one subgraph isomorphic to path structure in molecule
        # Store relative frequency of each structure (induced by path in rep_paths) as value in dictionary 'rep_strut'
        # Sort 'rep_strut' based on frequency highest to lowest
        # OrderedDict(sorted(rep_strut.items(), key=lambda x: x[1], reverse=True))
        return {}     # Return the sorted collection of path structures

    def add_structure_to_characteristic(self, sorted_list):
        # Use structures stored in order of frequency to decide where it should be added to characteristic substructure
        pass

    def find_characteristic_substructure(self, smiles_set):
        self.find_graphs_paths(smiles_set)
        length = self.length_start
        while length >= self.length_end:
            representative_paths = self.find_representative_paths(length)
            sorted_list = self.find_representative_structures(representative_paths)
            # After considering paths of this length test to see if there are representative substructures
            # If there are no rep structures then decrease stepwise, if there is increase the step size
            if sorted_list:
                self.add_structure_to_characteristic(sorted_list)
                length -= self.step
            else:
                length -= 1

if __name__ == '__main__':
    path_finder = CharacteristicSubstructure()
    print path_finder.find_graphs_paths(['CNO', 'POS', 'ClFIN', 'POSN'])
    print path_finder.find_representative_paths(3)
