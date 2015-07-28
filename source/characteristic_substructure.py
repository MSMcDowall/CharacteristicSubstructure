# coding=utf-8
from parser import Parser
import molecule
import subprocess
import networkx as nx
from draw_molecule import draw_molecule as draw


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
        self.paths = {}                     # All the paths from all the molecules with their lengths
        self.path_structures = {}           # Dictionary of path structures which has record of isomorphic molecule
        self.structure_text = {}

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
            if self.paths[path] == length:      # Check all the paths which are of the chosen length
                counter = 0
                for mole in self.molecules:
                    if path in (pair[0] for pair in mole.paths):        # Search tuples of (path, position) in molecule
                        counter += 1
                if float(counter)/len(self.molecules) >= self.path_threshold:
                    representative_paths.append(path)
        return representative_paths

    def create_structure(self, path, mole, vertices):
        structure = molecule.Molecule(path)
        molecule_map = {}
        for atom in vertices:
            if isinstance(atom, molecule.Atom):
                new_atom = structure.add_atom(atom.element)
                molecule_map[atom] = new_atom
            elif isinstance(atom, molecule.AromaticAtom):
                new_atom = structure.add_aromatic_atom(atom.element)
                molecule_map[atom] = new_atom
        for atom in vertices:
            for neighbour in vertices:
                edge = mole.contains_edge(atom, neighbour)      # Test if there is an edge to other atoms in path
                if isinstance(edge, molecule.SingleBond):
                    structure.add_single_bond(molecule_map[atom], molecule_map[neighbour])
                elif isinstance(edge, molecule.DoubleBond):
                    structure.add_double_bond(molecule_map[atom], molecule_map[neighbour])
                elif isinstance(edge, molecule.TripleBond):
                    structure.add_triple_bond(molecule_map[atom], molecule_map[neighbour])
                elif isinstance(edge, molecule.QuadrupleBond):
                    structure.add_quadruple_bond(molecule_map[atom], molecule_map[neighbour])
                elif isinstance(edge, molecule.AromaticBond):
                    structure.add_aromatic_bond(molecule_map[atom], molecule_map[neighbour])
        for vertex in structure.vertices:
            print vertex
        # Test if the structure has already been encountered
        self.create_text(structure)
        self.LAD_isomorphism(structure, mole, vertices)
        return structure

    def multiple_subgraphs(self, path, mole, positions):
        pass

    def create_text(self, path_structure):
        text = [str(path_structure.size) + '\n']
        for vertex in path_structure.vertices:
            for letter in vertex.element:
                number = ord(letter)
                text.append(str(number))
            text.append(' ' + str(len(path_structure.neighbours(vertex))) + ' ')
            for neighbour in path_structure.neighbours(vertex):
                text.append(str(neighbour.position) + ' ')
                edge = path_structure._vertices[vertex][neighbour]
                if isinstance(edge, molecule.SingleBond):
                    text.append('1 ')
                elif isinstance(edge, molecule.DoubleBond):
                    text.append('2 ')
                elif isinstance(edge, molecule.TripleBond):
                    text.append('3 ')
                elif isinstance(edge, molecule.QuadrupleBond):
                    text.append('4 ')
                elif isinstance(edge, molecule.AromaticBond):
                    text.append('5 ')
            text.append('\n')
        structure_string = ''.join(text)
        self.structure_text[path_structure] = structure_string
        print 'end of text'

    def LAD_isomorphism(self, pattern, mole, vertices):
        print 'entered isomorphism'
        pattern_file = open('pattern', mode='wb')
        pattern_file.write(self.structure_text[pattern])
        pattern_file.close()

        temporary_structure_dict = self.path_structures.copy()
        print self.path_structures.keys()
        print self.path_structures
        for structure in self.path_structures.keys():
            print' herey here'
            target_file = open('target', mode='wb')
            target_file.write(self.structure_text[structure])
            target_file.close()
            result = subprocess.check_output(["directedLAD/main", "-p", "pattern", "-t", "target", "-l", "-f"])
            result_location = 15    # Location in the string for the flag of success, 1 for success, 0 for failure
            print result
            if result[result_location] == '1':   # A match has been found and the mole and vertices have been added to that entry
                print type(self.path_structures[structure])
                temporary_structure_dict[structure][mole] = vertices
                print temporary_structure_dict
                break
        else:
            temporary_structure_dict[pattern] = {mole: vertices}
        self.path_structures = temporary_structure_dict.copy()
        print self.path_structures
        

    def find_representative_structures(self, rep_paths):
        for path in rep_paths:
            print path
            for mole in self.molecules:
                print mole.paths
                path_vertices = [pair[1] for pair in mole.paths if pair[0] == path]    # Lists of path vertices
                print 'vert in path'
                print path_vertices
                # If path_vertices > 1 go to special method
                if len(path_vertices) > 0:
                    if len(path_vertices) > 1:
                        self.multiple_subgraphs(path, mole, path_vertices)
                    else:
                        print 'len greater than 0'
                        self.create_structure(path, mole, path_vertices[0])
        # Consider here what happens if there is more than one subgraph isomorphic to path structure in molecule
        # Store relative frequency of each structure (induced by path in rep_paths) as value in dictionary 'rep_strut'
        # Sort 'rep_strut' based on frequency highest to lowest
        # OrderedDict(sorted(rep_strut.items(), key=lambda x: x[1], reverse=True))

    # NOT YET IN USE
    def add_structure_to_characteristic(self, sorted_list):
        # Use structures stored in order of frequency to decide where it should be added to characteristic substructure
        pass

    # NOT YET IN USE
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
    path_finder = CharacteristicSubstructure(path_threshold=0.3)
    print 'all paths'
    print path_finder.find_graphs_paths(['CNO', 'ONC', 'CNO'])
    rep_paths = path_finder.find_representative_paths(3)
    print 'rep paths'
    print rep_paths
    path_finder.find_representative_structures(rep_paths)
    print path_finder.path_structures
    print 'keys!'
    for strut in path_finder.path_structures:
        print strut
        #print path_finder.path_structures[strut]
    #draw(path_finder.path_structures.keys()[0])
    #draw(path_finder.path_structures.keys()[1])

    # tester = CharacteristicSubstructure()
    # molecule1 = Parser().parse_smiles('CNO')
    # # molecule2 = Parser().parse_smiles('ONC')
    # # molecule3 = Parser().parse_smiles('FPC')
    # tester.create_text(molecule1)
    # # tester.create_text(molecule2)
    # # tester.create_text(molecule3)
    # tester.path_structures[molecule1] = {'first': 'initial value1'}
    # # tester.path_structures[molecule2] = {'second': 'initial value 2'}
    # # tester.path_structures[molecule3] = {'third': 'initial value 3'}
    # tester.LAD_isomorphism(molecule1, 'first', ['new value'])
    # # tester.LAD_isomorphism(molecule2, 'second', ['2nd value'])
    # # tester.LAD_isomorphism(molecule3, 'third', ['3rd value'])


