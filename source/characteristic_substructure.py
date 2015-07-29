# coding=utf-8
from parser import Parser
import molecule
import subprocess
from collections import OrderedDict
import networkx as nx
import networkx.algorithms.isomorphism as iso
from draw_molecule import draw_molecule as draw


# Implementation of Finding Characteristic Substructures for Metabolite Classes
# Ludwig, Hufsky, Elshamy, Böcker


class CharacteristicSubstructure(object):
    def __init__(self, length_start=20, length_end=5, step=4, threshold=0.8, isomorphism_factor=0.8):
        self.length_start = length_start
        self.length_end = length_end
        self.step = step
        self.threshold = threshold
        self.isomorphism_factor = isomorphism_factor
        self.molecules = []                 # All the given molecules
        self.paths = {}                     # All the paths from all the molecules with their lengths
        self.path_structures = {}           # Dictionary of path structures which has record of isomorphic molecule
        self.structure_text = {}            # Text version of graph for LAD isomorphism
        self.structure_nx = {}              # NetworkX graph version of structure
        self.multiple_structures = {}       # Dictionary of the structures which appear multiple times in molecules

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
                if float(counter)/len(self.molecules) >= self.threshold:
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
        return structure

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

    # ERROR Treats rings as directed graph
    def lad_isomorphism(self, pattern, mole, vertices):
        print 'entered LAD isomorphism'
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

    def create_nx_graph(self, path_structure):
        # Create a new NetworkX graph
        g = nx.Graph()
        # For each vertex and edge in molecule graph add node and edge in NetworkX graph
        for n in path_structure.vertices():
            g.add_node(n.position, element=n.element)
        for e in path_structure.edges():
            if isinstance(e, molecule.SingleBond):
                g.add_edge(e.endpoints_position()[0], e.endpoints_position()[1], type='single')
            elif isinstance(e, molecule.DoubleBond):
                g.add_edge(e.endpoints_position()[0], e.endpoints_position()[1], type='double')
            elif isinstance(e, molecule.TripleBond):
                g.add_edge(e.endpoints_position()[0], e.endpoints_position()[1], type='triple')
            elif isinstance(e, molecule.QuadrupleBond):
                g.add_edge(e.endpoints_position()[0], e.endpoints_position()[1], type='quadruple')
            elif isinstance(e, molecule.AromaticBond):
                g.add_edge(e.endpoints_position()[0], e.endpoints_position()[1], type='aromatic')
        self.structure_nx[path_structure] = g

    def nx_isomorphism(self, pattern, mole, vertices):
        print 'entered nx isomorphism'
        temporary_structure_dict = self.path_structures.copy()
        for structure in self.path_structures.keys():
            if not nx.faster_could_be_isomorphic(self.structure_nx[pattern], self.structure_nx[structure]):
                # Graphs are definitely not isomorphic
                continue    # Continues to next structure for testing
            matcher = iso.GraphMatcher(self.structure_nx[pattern], self.structure_nx[structure],
                                       node_match=iso.categorical_node_match('element', 'C'),
                                       edge_match=iso.categorical_edge_match('type', 'single'))
            if matcher.is_isomorphic():
                temporary_structure_dict[structure][mole] = vertices
                # if mole in temporary_structure_dict[structure]:
                #     # There are multiple occurrences of structure in molecule
                #     if structure in self.multiple_structures and mole in self.multiple_structures[structure]:
                #         self.multiple_structures[structure][mole].append(vertices)
                #     else:
                #         self.multiple_structures[structure] = {mole: [vertices]}
                #     # As the structure occurs multiple times we do not want to choose
                #     del temporary_structure_dict[structure][mole]
                # else:
                #     temporary_structure_dict[structure] = {mole: vertices}
                break
        else:
            temporary_structure_dict[pattern] = {mole: vertices}
        self.path_structures = temporary_structure_dict.copy()
        print self.path_structures

    def find_representative_structures(self, rep_paths):
        rep_structures = {}
        for path in rep_paths:
            print path
            for mole in self.molecules:
                print mole.paths
                path_vertices = [pair[1] for pair in mole.paths if pair[0] == path]    # Lists of path vertices
                print 'vert in path'
                print path_vertices
                for vertices in path_vertices:
                    structure = self.create_structure(path, mole, vertices)
                    # Test if the structure has already been encountered - change here between LAD and nx isomorphism
                    # self.create_text(structure)
                    # self.LAD_isomorphism(structure, mole, vertices)
                    self.create_nx_graph(structure)
                    self.nx_isomorphism(structure, mole, vertices)
        for structure in self.path_structures:
            print 'length'
            print len(self.path_structures[structure].keys())
            print len(self.molecules)
            relative_frequency = len(self.path_structures[structure].keys())/float(len(self.molecules))
            print float(relative_frequency)
            if float(relative_frequency) >= self.threshold:
                rep_structures[structure] = relative_frequency
        # Store relative frequency of each structure (induced by path in rep_paths) as value in dictionary
        # Sort dictionary based on frequency highest to lowest
        representative_structures = OrderedDict(sorted(rep_structures.items(), key=lambda x: x[1], reverse=True))
        return representative_structures

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
    path_finder = CharacteristicSubstructure(threshold=0.3)
    print 'all paths'
    print path_finder.find_graphs_paths(['CNO', 'CNOP', 'FSP'])
    rep_paths = path_finder.find_representative_paths(2)
    print 'rep paths'
    print rep_paths
    rep_struct = path_finder.find_representative_structures(rep_paths)
    print 'rep struct'
    print rep_struct



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


