# coding=utf-8
from smiles_parser import Parser
import molecule as m
from collections import OrderedDict, Counter
from copy import copy
import networkx as nx
import networkx.algorithms.isomorphism as iso
# from draw_molecule import draw_molecule as draw


# Implementation of Finding Characteristic Substructures for Metabolite Classes
# Ludwig, Hufsky, Elshamy, Böcker

class CharacteristicSubstructure(object):
    def __init__(self, length_start=20, length_end=5, step=4, threshold=0.8, isomorphism_factor=0.8):
        # The initial parameters for the algorithm
        self.length_start = length_start
        self.length_end = length_end
        self.step = step
        self.threshold = threshold
        self.isomorphism_factor = isomorphism_factor
        self.characteristic_substructure = None
        self.molecules = []                 # All the given molecules
        self.paths = {}                     # All the paths from all the molecules with their lengths
        self.path_structures = {}           # Dictionary of path structures which has record of isomorphic molecules
        self.structure_nx = {}              # NetworkX graph version of structure
        self.multiple_structures = {}       # Dictionary of the paths which appear multiple times in molecules

    def find_graphs_paths(self, smiles_set):
        for smiles in smiles_set:
            mole = Parser().parse_smiles(smiles)
            self.molecules.append(mole)
            path_dict = mole.find_all_paths()
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

    def find_representative_structures(self, rep_paths):
        rep_structures = {}
        for path in rep_paths:
            # print path
            for mole in self.molecules:
                path_vertices = [pair[1] for pair in mole.paths if pair[0] == path]    # Lists of path vertices
                for vertices in path_vertices:
                    structure_tuple = self._create_structure(path, mole, vertices)
                    structure = structure_tuple[0]
                    vertices = structure_tuple[1]
                    # Test if the structure has already been encountered - change here between LAD and nx isomorphism
                    self._create_nx_graph(structure)
                    self._check_structure_duplicates(structure, mole, vertices)
        # print self.path_structures
        if self.multiple_structures:
            self._create_multiple_structures()
        for structure in self.path_structures:
            relative_frequency = len(self.path_structures[structure].keys())/float(len(self.molecules))
            if float(relative_frequency) >= self.threshold:
                rep_structures[structure] = relative_frequency
        # Store relative frequency of each structure (induced by path in rep_paths) as value in dictionary
        # Sort dictionary based on frequency highest to lowest
        representative_structures = OrderedDict(sorted(rep_structures.items(), key=lambda x: x[1], reverse=True))
        return representative_structures

    def find_characteristic_substructure(self):
        smiles_set = self._data_input()
        self.find_graphs_paths(smiles_set)
        length = self.length_start
        while length >= self.length_end:
            representative_paths = self.find_representative_paths(length)
            sorted_list = self.find_representative_structures(representative_paths)
            # After considering paths of this length test to see if there are representative substructures
            # If there are no rep structures then decrease stepwise, if there is increase the step size
            if sorted_list:
                print 'sorted list of representative structures'
                print sorted_list
                for structure in sorted_list:
                    # draw(structure)
                    self._add_structure_to_characteristic(structure)
                length -= self.step
            else:
                length -= 1
        return self.characteristic_substructure

    def _add_structure_to_characteristic(self, structure):
        # When a structure is added to the CS, any molecule that is sub isomorphic to it will have its vertices swapped
        # These molecules will instead contain CS vertices
        print 'structure to be added to CS'
        # print self.path_structures[structure]
        print structure.adjacency_dictionary.keys()
        if not self.characteristic_substructure:
            self.characteristic_substructure = copy(structure)
            print 'CS made'
            print self.characteristic_substructure.adjacency_dictionary.keys()
            # draw(self.characteristic_substructure)
        elif self.characteristic_substructure:
            print 'adding to CS'
            print str(structure)
            possible_locations = []    # List of possible locations of subgraphs
            for mole in self.path_structures[structure]:
                # Create copy of the structure so that changes can be made without altering original structure
                structure_copy = m.Molecule(str(structure))
                structure_copy.adjacency_dictionary = structure.adjacency_dictionary.copy()
                # Does structure vertices map to any CS vertices?
                # If they do then change them to the CS vertices
                for vertex in self.path_structures[structure][mole]:
                    print 'switching out copy of structure'
                    print repr(vertex)
                    print structure_copy.adjacency_dictionary[vertex]
                    print structure_copy.adjacency_dictionary.keys()
                    if self.path_structures[structure][mole][vertex] in self.characteristic_substructure.adjacency_dictionary:
                        structure_copy.swap_vertex(vertex, self.path_structures[structure][mole][vertex])
                    print 'successful switch'
                # print 'new altered struct vertices'
                # print structure.adjacency_dictionary.keys()
                # print structure_copy.adjacency_dictionary.keys()
                # Create a copy of the characteristic substructure which will have the structure added to it
                possible_location = m.Molecule(str(self.characteristic_substructure))
                possible_location.adjacency_dictionary = self.characteristic_substructure.adjacency_dictionary.copy()
                # print 'cs vertices'
                # print self.characteristic_substructure.adjacency_dictionary.keys()
                # print possible_location.adjacency_dictionary.keys()
                # For vertices in structure not in CS
                # Add to CS also append bonds for structure vertices
                for vertex in structure_copy.adjacency_dictionary:
                    if vertex in self.characteristic_substructure.adjacency_dictionary:
                        possible_location.adjacency_dictionary[vertex].update(structure_copy.adjacency_dictionary[vertex])
                    else:
                        possible_location.adjacency_dictionary[vertex] = structure_copy.adjacency_dictionary[vertex]
                print 'added to adj dictionary'
                # print possible_location.adjacency_dictionary
                # draw(possible_location)
                possible_locations.append(possible_location)
            print possible_locations
            for possible in possible_locations:
                self._create_nx_graph(possible)
            isomorphic_locations = {}
            for possible in possible_locations:
                # print 'checking locations for iso'
                # print isomorphic_locations
                for location in isomorphic_locations:
                    if self._nx_isomorphism(possible, location):
                        isomorphic_locations[location] += 1
                        break
                else:
                    isomorphic_locations[possible] = 1
            print isomorphic_locations
            sorted_locations = OrderedDict(sorted(isomorphic_locations.items(), key=lambda x: x[1], reverse=True))
            self.characteristic_substructure = sorted_locations.items()[0][0]
            print 'new cs'
            print self.characteristic_substructure
            # draw(self.characteristic_substructure)
        print 'CS vertices'
        print self.characteristic_substructure.adjacency_dictionary.keys()
        self._swap_molecule_vertices(structure)

    def _create_structure(self, path, mole, vertices):
        structure = m.Molecule(path)
        old_molecule_map = {}           # Uses the original vertices as keys, structure vertices as values
        new_molecule_map = {}           # Uses the structure vertices as keys, original vertices as values
        for atom in vertices:
            if isinstance(atom, m.Atom):
                new_atom = structure.add_atom(atom.element)
                old_molecule_map[atom] = new_atom
                new_molecule_map[new_atom] = atom
            elif isinstance(atom, m.AromaticAtom):
                new_atom = structure.add_aromatic_atom(atom.element)
                old_molecule_map[atom] = new_atom
                new_molecule_map[new_atom] = atom
        for atom in vertices:
            for neighbour in vertices:
                # ERROR atom key error
                edge = mole.contains_edge(atom, neighbour)      # Test if there is an edge to other atoms in path
                if isinstance(edge, m.SingleBond):
                    structure.add_single_bond(old_molecule_map[atom], old_molecule_map[neighbour])
                elif isinstance(edge, m.DoubleBond):
                    structure.add_double_bond(old_molecule_map[atom], old_molecule_map[neighbour])
                elif isinstance(edge, m.TripleBond):
                    structure.add_triple_bond(old_molecule_map[atom], old_molecule_map[neighbour])
                elif isinstance(edge, m.QuadrupleBond):
                    structure.add_quadruple_bond(old_molecule_map[atom], old_molecule_map[neighbour])
                elif isinstance(edge, m.AromaticBond):
                    structure.add_aromatic_bond(old_molecule_map[atom], old_molecule_map[neighbour])
        return structure, new_molecule_map

    def _create_multiple_structures(self):
        # multiple_structures = {structure: mole: [[strut, vert map]...]}
        for structure in self.multiple_structures:
            # "If at least k>= 1 isomorphic subgraphs within the same molecule graph in at least |M|.iso graphs occur:"
            # If the subgraph structure appears in |M|.iso graphs then multiple structures will be created
            if len(self.path_structures[structure]) > len(self.molecules)*self.isomorphism_factor:
                for mole in self.multiple_structures[structure]:
                    k = len(self.multiple_structures[structure][mole])
                    path = (str(structure) + ' ') * k
                    new_structure = m.Molecule(path)
                    new_structure.adjacency_dictionary.update(structure.adjacency_dictionary)
                    new_structure.size += structure.size
                    vertices_mapping = self.path_structures[structure][mole].copy()
                    for pair in self.multiple_structures[structure][mole]:
                        substructure = pair[0]
                        # Ensure positions in the structure are unique so they can be used with the nx graph
                        for vertex in substructure.adjacency_dictionary:
                            vertex.position += new_structure.size
                        new_structure.adjacency_dictionary.update(substructure.adjacency_dictionary)
                        new_structure.size += substructure.size
                        vertices_mapping.update(pair[1])
                    # Remove that mole entry from path structures dictionary so it is not repeated
                    del self.path_structures[structure][mole]
                    # Add the new larger structure into the dictionary
                    self._create_nx_graph(new_structure)
                    self._check_structure_duplicates(new_structure, mole, vertices_mapping)
                    # draw(new_structure)

    def _create_nx_graph(self, path_structure):
        # Create a new NetworkX graph
        g = nx.Graph()
        # For each vertex and edge in molecule graph add node and edge in NetworkX graph
        for n in path_structure.vertices():
            g.add_node(n.position, element=n.element)
        for e in path_structure.edges():
            if isinstance(e, m.SingleBond):
                g.add_edge(e.endpoints_position()[0], e.endpoints_position()[1], type='single')
            elif isinstance(e, m.DoubleBond):
                g.add_edge(e.endpoints_position()[0], e.endpoints_position()[1], type='double')
            elif isinstance(e, m.TripleBond):
                g.add_edge(e.endpoints_position()[0], e.endpoints_position()[1], type='triple')
            elif isinstance(e, m.QuadrupleBond):
                g.add_edge(e.endpoints_position()[0], e.endpoints_position()[1], type='quadruple')
            elif isinstance(e, m.AromaticBond):
                g.add_edge(e.endpoints_position()[0], e.endpoints_position()[1], type='aromatic')
        self.structure_nx[path_structure] = g

    def _nx_isomorphism(self,pattern, target):
        if not nx.faster_could_be_isomorphic(pattern, target):
                # Graphs are definitely not isomorphic
                return None
        matcher = iso.GraphMatcher(pattern, target,
                                   node_match=iso.categorical_node_match('element', 'C'),
                                   edge_match=iso.categorical_edge_match('type', 'single'))
        if matcher.is_isomorphic():
            return matcher.mapping

    def _check_structure_duplicates(self, pattern, mole, vertices):
        temporary_structure_dict = self.path_structures.copy()
        for structure in self.path_structures.keys():
            isomorphic_mapping = self._nx_isomorphism(self.structure_nx[pattern], self.structure_nx[structure])
            if not isomorphic_mapping:
                continue    # Continues to next structure for testing as they are not isomorphic
            if isomorphic_mapping:
                if mole in temporary_structure_dict[structure]:
                    # If the molecule vertices are different to the vertices currently in dictionary
                    # then there are multiple subgraphs in molecule isomorphic to pattern
                    if Counter(vertices.values()) != Counter(temporary_structure_dict[structure][mole].values()):
                        if structure in self.multiple_structures and mole in self.multiple_structures[structure]:
                            self.multiple_structures[structure][mole].append([pattern, vertices])
                        else:
                            self.multiple_structures[structure] = {mole: [[pattern, vertices]]}
                else:
                    altered_mapping = {}
                    for g2_position in isomorphic_mapping:
                        for vertex in structure.adjacency_dictionary:
                            if vertex.position == g2_position:
                                g2_match = vertex
                                break
                        for vertex in pattern.adjacency_dictionary:
                            if vertex.position == isomorphic_mapping[g2_position]:
                                mole_vertex = vertices[vertex]
                                break
                        altered_mapping[g2_match] = mole_vertex
                    temporary_structure_dict[structure][mole] = altered_mapping
                break
        else:
            temporary_structure_dict[pattern] = {mole: vertices}
        self.path_structures = temporary_structure_dict.copy()

    def _data_input(self):
        smiles_set = []
        reader = open('SMILES.txt', mode='rb')
        for line in reader:
            smiles_set.append(line.rstrip())
        reader.close()
        return smiles_set

    def _data_output(self):
        pass

    def _swap_molecule_vertices(self, structure):
        for mole in self.path_structures[structure]:
            for key in self.path_structures[structure][mole]:
                # If this vertex of the structure has been added to the CS
                # then change the molecule that are associated with it
                if key in self.characteristic_substructure.adjacency_dictionary.keys():
                    # Change the vertices in the molecule
                    # ERROR old_vertex invalid key
                    mole.swap_vertex(self.path_structures[structure][mole][key], key)
                    self._swap_path_structure(self.path_structures[structure][mole][key], key)

    def _swap_path_structure(self, old, new):
        for structure in self.path_structures:
            for molecule in self.path_structures[structure]:
                for vertex in self.path_structures[structure][molecule]:
                    if self.path_structures[structure][molecule][vertex] == old:
                        self.path_structures[structure][molecule][vertex] = new

if __name__ == '__main__':
    path_finder = CharacteristicSubstructure(threshold=0.3)
    structure = path_finder.find_characteristic_substructure()
    print structure
