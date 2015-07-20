from parser import Parser
# This will contain the methods to find paths in all molecules
# Test if path is in graph
# Check frequency of each path


class RepresentativePaths(object):
    def __init__(self):
        self.paths = {}
        self.molecules = []
        self.path_threshold = 0.5
        self.representative_paths = []

    def find_graphs_paths(self, smiles_set):
        for smiles in smiles_set:
            molecule = Parser().parse_smiles(smiles)
            self.molecules.append(molecule)
            path_dict = molecule.find_all_paths()
            self.paths.update(path_dict)
        return self.paths

    def find_representative_paths(self, length):
        for path in self.paths:
            if self.paths[path] == length:
                counter = 0
                for mole in self.molecules:
                    if path in mole.paths:
                        counter += 1
                if float(counter)/len(self.molecules) >= self.path_threshold:
                    self.representative_paths.append(path)
        return self.representative_paths


if __name__ == '__main__':
    path_finder = RepresentativePaths()
    print path_finder.find_graphs_paths(['CNO', 'POS', 'ClFIN', 'POSN'])
    print path_finder.find_representative_paths(3)