# coding=utf-8
from smiles_parser import Parser
import molecule as m
from collections import OrderedDict, Counter
from copy import copy
import sys
import networkx as nx
import networkx.algorithms.isomorphism as iso

# Implementation of Finding Characteristic Substructures for Metabolite Classes
# Ludwig, Hufsky, Elshamy, Böcker


class CSAlgorithm(object):
    def __init__(self, smiles_file="SMILES.txt", length_start=20, length_end=5, threshold=0.8):
        # The file containing the SMILES strings
        self.smiles_file = smiles_file
        # The initial parameters for the algorithm
        self.length_start = length_start
        self.length_end = length_end
        self.step = 4
        self.threshold = threshold
        self.isomorphism_factor = 0.8

        # The structure which has been created through the combination of representative path structures
        self.characteristic_substructure = m.Molecule('Characteristic Substructure')
        # List holding structures which have been added to the characteristic substructure
        self.cs_structures = []
        # Dictionary holding the locations of the molecules which map to the characteristic substructure
        self.cs_locations = {}
        # All the given molecules
        self.molecules = []
        # Dictionary of path structures with the isomorphic molecules and mapping vertices from structure to molecule
        # {path structure: {molecule: {structure vertex: molecule vertex}}}
        self.path_structures = {}
        # Dictionary of structures that appear multiple times in molecules and list of molecule vertices that map them
        # {single structure: {molecule: {structure vertex: [molecule vertices]}}}
        self.multiple_vertices = {}
        # Dictionary of structures which are composed of instances of other path structures, maps to molecule vertices
        # {multiple structure: {molecule: {structure vertex: molecule vertex}}}
        self.multiple_structures = {}
        # NetworkX graph version of the path structure, path structure is key
        self.structure_nx = {}

    def find_characteristic_substructure(self):
        """
        Find the characteristic substructure for a set of molecules

        Calls the representative paths method for each of the lengths between the start length and the end
        Creates the CS with the most frequent path structure
        Each of the subsequent path structures is added to the CS in order of frequency
        Swaps the vertices in the molecules which map to the CS with the CS vertices

        :return: a molecule object that is the characteristic substructure of the list of molecules
        """
        smiles_set = self._data_input()
        paths = self._find_graphs_paths(smiles_set)
        length = self.length_start
        while length >= self.length_end:
            print length
            representative_paths = self._find_representative_paths(paths, length)
            print representative_paths
            sorted_dictionary = self._find_representative_structures(representative_paths)
            print sorted_dictionary
            # After considering paths of this length test to see if there are representative substructures
            # If there are no rep structures then decrease stepwise, if there is increase the step size
            if sorted_dictionary:
                for structure in sorted_dictionary.keys():
                    print 'adding...'
                    self._add_structure_to_characteristic(structure)
                length -= self.step
            else:
                length -= 1
        self._data_output(self._create_cs_results(), self.smiles_file[:-4] + str(self.threshold) +
                          '_CharacteristicSubstructure.txt')
        return self.characteristic_substructure

    def find_all_representative_structures(self):
        """
        Creates a list of all the structure of different lengths which are representative sorted in terms of frequency

        :return: list of structures which appear frequently in molecules
        """
        all_structures = {}
        smiles_set = self._data_input()
        paths = self._find_graphs_paths(smiles_set)
        length = self.length_start
        while length >= self.length_end:
            representative_paths = self._find_representative_paths(paths, length)
            sorted_dictionary = self._find_representative_structures(representative_paths)
            # After considering paths of this length test to see if there are representative substructures
            # If there are no rep structures then decrease stepwise, if there is increase the step size
            if sorted_dictionary:
                all_structures.update(sorted_dictionary)
            # To get the structures of all lengths the step does not alter
            length -= 1
        representative_structures = OrderedDict(sorted(all_structures.items(), key=lambda x: x[1], reverse=True)).keys()
        self._data_output(self._structures_output(representative_structures),
                          self.smiles_file[:-4] + str(self.threshold) + '_RepresentativeStructures.txt')
        return representative_structures

    def _find_graphs_paths(self, smiles_set):
        """
        For each SMILES string a molecule object is created and all of its paths found
        
        :param smiles_set: a list of SMILES strings
        :return: a dictionary containing all of the paths as strings with their length as value
        """
        paths = {}
        for smiles in smiles_set:
            molecule = Parser().parse_smiles(smiles)
            # draw(molecule)
            self.molecules.append(molecule)
            path_dict = molecule.find_all_paths()
            paths.update(path_dict)
        return paths

    def _find_representative_paths(self, paths, length):
        """
        Find the paths which occur with a high enough frequency
        
        :param length: an integer which sets the length of the paths that should be considered
        :return: a list of the paths as strings which are representative 
        """
        representative_paths = []
        for path in paths:
            # Check all the paths which are of the chosen length
            if paths[path] == length:
                counter = 0
                for molecule in self.molecules:
                    # Search the tuples for each path which consists of path string and vertices present in path
                    if path in (pair[0] for pair in molecule.paths):
                        counter += 1
                if float(counter) / len(self.molecules) >= self.threshold:
                    representative_paths.append(path)
        return representative_paths

    def _find_representative_structures(self, rep_paths):
        """
        Find the path structures which appear with a high enough frequency.

        Ensures there are no duplicates in the dictionary of path structures so that the frequency count is accurate
        :param rep_paths: a list of the paths which are representative
        :return: a list of path structures which are representative sorted into order of frequency, highest to lowest
        """
        rep_structures = {}
        for path in rep_paths:
            for molecule in self.molecules:
                path_vertices = [pair[1] for pair in molecule.paths if pair[0] == path]  # Lists of path vertices
                for vertices in path_vertices:
                    structure_tuple = self._create_structure(path, molecule, vertices)
                    structure = structure_tuple[0]
                    vertices = structure_tuple[1]
                    # Test if the structure has already been encountered
                    self._check_structure_duplicates(structure, molecule, vertices)
        for structure in self.path_structures:
            relative_frequency = len(self.path_structures[structure].keys()) / float(len(self.molecules))
            if float(relative_frequency) >= self.threshold:
                # Store relative frequency of each structure (induced by path in rep_paths) as value in dictionary
                rep_structures[structure] = relative_frequency
        # Sort dictionary based on frequency highest to lowest and return the path structures only
        return OrderedDict(sorted(rep_structures.items(), key=lambda x: x[1], reverse=True))

    def _create_structure(self, path, molecule, vertices):
        """
        Creates a path structure which is a molecule object representing the form that the given path takes within a molecule

        :param path: the string which represents the path (Entered as the SMILES string for the path structure)
        :param molecule: the molecule object which contains the given path and vertices
        :param vertices: a list of vertices which are encountered in the path of the molecule
        :return: a tuple of the newly created path structure and
                a dictionary mapping the structures vertices to the molecule vertices
        """
        structure = m.Molecule(path)
        # Uses the original vertices as keys, structure vertices as values, used to create the edges
        old_molecule_map = {}
        # Uses structure vertices as keys, original vertices as values,
        # used to find molecule vertices if structure is given
        new_molecule_map = {}
        for atom in vertices:
            if not atom.aromatic:
                new_atom = structure.add_atom(atom.element)
                old_molecule_map[atom] = new_atom
                new_molecule_map[new_atom] = atom
            elif atom.aromatic:
                new_atom = structure.add_aromatic_atom(atom.element)
                old_molecule_map[atom] = new_atom
                new_molecule_map[new_atom] = atom
        for atom in vertices:
            for neighbour in vertices:
                # Test if there is an edge to any of the other atoms in the path
                edge = molecule.contains_edge(atom, neighbour)
                if isinstance(edge, m.Bond):
                    if edge.single:
                        structure.add_single_bond(old_molecule_map[atom], old_molecule_map[neighbour])
                    elif edge.double:
                        structure.add_double_bond(old_molecule_map[atom], old_molecule_map[neighbour])
                    elif edge.triple:
                        structure.add_triple_bond(old_molecule_map[atom], old_molecule_map[neighbour])
                    elif edge.quadruple:
                        structure.add_quadruple_bond(old_molecule_map[atom], old_molecule_map[neighbour])
                    elif edge.aromatic:
                        structure.add_aromatic_bond(old_molecule_map[atom], old_molecule_map[neighbour])
        return structure, new_molecule_map

    def _check_structure_duplicates(self, pattern, molecule, vertices):
        """
        Adds the pattern graph to path_structures or multiple_structures depending on if there is a duplicate of it.

        :param pattern: the graph object which is to be tested against the other graphs
        :param molecule: the molecule which contained the subgraph isomorphic to the structure
        :param vertices: the dictionary mapping the pattern vertices to the molecule vertices
        :return: the graph object which has been newly added or updated
        """
        unique_structure = None
        # A copy of path_structures is made so that it can be altered while also looping over the element
        temporary_structure_dict = self.path_structures.copy()
        for structure in self.path_structures:
            nx_mapping = self._nx_isomorphism(pattern, structure)
            if not nx_mapping:
                # Continues to next structure for testing as they are not isomorphic
                continue
            if nx_mapping:
                # A multiple structure has been made which matches a structure that has already been added to the CS
                if structure in self.cs_structures:
                    return False
                # Maps the vertices from the pattern to the target
                isomorphic_mapping = {}
                for target_position in nx_mapping:
                    for vertex in structure.adjacency_dictionary:
                        if vertex.position == target_position:
                            target_match = vertex
                            break
                    for vertex in pattern.adjacency_dictionary:
                        if vertex.position == nx_mapping[target_position]:
                            mole_vertex = vertices[vertex]
                            break
                    isomorphic_mapping[target_match] = mole_vertex
                if molecule in temporary_structure_dict[structure]:
                    # If the molecule vertices are different to the vertices currently in dictionary
                    # then there are multiple subgraphs in molecule isomorphic to pattern
                    if Counter(isomorphic_mapping.values()) != Counter(temporary_structure_dict[structure][molecule].values()):
                        self._add_structure_to_multiple_dictionary(structure, molecule, isomorphic_mapping)
                else:
                    temporary_structure_dict[structure][molecule] = isomorphic_mapping
                unique_structure = structure
                break
        else:
            temporary_structure_dict[pattern] = {molecule: vertices}
            unique_structure = pattern
        self.path_structures = temporary_structure_dict.copy()
        return unique_structure

    def _nx_isomorphism(self, pattern, target):
        """
        Uses the NetworkX isomorphism algorithm to check if the pattern graph and the target graph are isomorphic.

        The faster_could_be_isomorphic method is used to discount two structures if they could not be isomorphic.
        :param pattern: a molecule object which is to be tested for isomorphism
        :param target: a molecule object which pattern graph is to be compared against
        :return: None if the graphs are not isomorphic
        :return: a dictionary which maps the indices of the two NetworkX graphs together if they are isomorphic
        """
        if pattern not in self.structure_nx:
            self._create_nx_graph(pattern)
        if target not in self.structure_nx:
            self._create_nx_graph(target)
        if not nx.faster_could_be_isomorphic(self.structure_nx[pattern], self.structure_nx[target]):
            # Graphs are definitely not isomorphic
            return None
        # Ensures the isomorphism considers the vertex element and edge type
        matcher = iso.GraphMatcher(self.structure_nx[pattern], self.structure_nx[target],
                                   node_match=iso.categorical_node_match('element', 'C'),
                                   edge_match=iso.categorical_edge_match('type', 'single'))
        if matcher.is_isomorphic():
            return matcher.mapping

    def _create_nx_graph(self, path_structure):
        """
        Creates a copy of the molecule object as an NetworkX graph
        The position attirbute of each structure vertices is used as the index of the NetworkX graph

        :param path_structure: the molecule object which will be turned into a NetworkX graph
        :return: None
        """
        g = nx.Graph()
        # For each vertex and edge in molecule graph add node and edge in NetworkX graph
        for n in path_structure.vertices():
            g.add_node(n.position, element=n.element)
        for e in path_structure.edges():
            if isinstance(e, m.Bond):
                if e.single:
                    g.add_edge(e.endpoints_position()[0], e.endpoints_position()[1], type='single')
                elif e.double:
                    g.add_edge(e.endpoints_position()[0], e.endpoints_position()[1], type='double')
                elif e.triple:
                    g.add_edge(e.endpoints_position()[0], e.endpoints_position()[1], type='triple')
                elif e.quadruple:
                    g.add_edge(e.endpoints_position()[0], e.endpoints_position()[1], type='quadruple')
                elif e.aromatic:
                    g.add_edge(e.endpoints_position()[0], e.endpoints_position()[1], type='aromatic')
        self.structure_nx[path_structure] = g

    def _add_structure_to_multiple_dictionary(self, structure, molecule, mapping):
        """
        Adds a structure which has appeared more than once in a molecule to the multiple_structures dictionary

        :param structure: the structure which has appeared more than once in a molecule
        :param molecule: the molecule that has repeated instances of the structure
        :param mapping: a mapping from the structure to the molecule vertices which comes from the isomorphism check
        :return: None
        """
        if structure in self.multiple_vertices and molecule in self.multiple_vertices[structure]:
            # Creates a list of lists which sets out the vertices that are mapped to each instance of the structure
            # within the molecule, this ensures that an instance isn't duplicated within the dictionary.
            multiples = [[vertex[index] for vertex in self.multiple_vertices[structure][molecule].values()]
                         for index in range(len(self.multiple_vertices[structure][molecule].values()[0]))]
            for v in multiples:
                if Counter(mapping.values()) == Counter(v):
                    return
            # Store the molecule vertex which maps to each of the structure vertices
            for key in mapping:
                self.multiple_vertices[structure][molecule][key].append(mapping[key])
        else:
            if structure in self.multiple_vertices:
                self.multiple_vertices[structure][molecule] = {}
            elif structure not in self.multiple_vertices:
                self.multiple_vertices[structure] = {molecule: {}}
            # When the structure is encountered for the second time in the molecule the mapping from path_structures
            # is also stored into the multiple_structures dictionary so it includes all instances of the structure
            for key in mapping:
                self.multiple_vertices[structure][molecule][key] = [self.path_structures[structure][molecule][key],
                                                                      mapping[key]]

    def _add_structure_to_characteristic(self, structure):
        """
        Adds the given structure to the characteristic substructure in the location where it appears most frequently
        A different method is called depending on if the structure appears once or many times in the molecules
        A list of the possible locations on the CS where the structure could be added is created
        This list is then sorted in terms of frequency and the CS with the most common position becomes the new CS

        :param structure: the molecule object which is to be added to the characteristic substructure
        :return: None
        """
        if structure in self.multiple_vertices:
            self._add_multiple_to_characteristic(structure)
        elif structure not in self.multiple_vertices:
            possible_locations = []
            for molecule in self.path_structures[structure]:
                possible = self._add_single_to_characteristic(structure, molecule)
                possible_locations.append(possible)
            self.characteristic_substructure = self._most_frequent_location(possible_locations)
            self._add_cs_locations(structure)

    def _add_single_to_characteristic(self, structure, molecule):
        """
        Creates one of the possible locations where the structure could be added to characteristic substructure.

        A copy of the structure is added to a graph which includes a copy of the CS
        These copies are used so that the original adjacency dictionaries of the structure and CS are not altered.
        :param structure: The structure that is to be added to CS in different locations
        :return: list of graphs which display the different locations where the structure could be added to CS
        """
        # Create a copy of the characteristic substructure which will have the structure added to it
        possible_location = m.Molecule(str(self.characteristic_substructure))
        for vertex in self.characteristic_substructure.adjacency_dictionary:
            possible_location.vertex_to_graph(vertex)
            for neighbour in self.characteristic_substructure.adjacency_dictionary[vertex]:
                possible_location.adjacency_dictionary[vertex][neighbour] = copy(self.characteristic_substructure.
                                                                                 adjacency_dictionary
                                                                                 [vertex][neighbour])
        # This tests if the this molecule has any subgraphs isomorphic to the current characteristic substructure
        # If there are no isomorphic subgraphs then the possible CS will be disconnected
        if molecule not in self.cs_locations:
            for vertex in structure.adjacency_dictionary:
                possible_location.vertex_to_graph(vertex)
                possible_location.adjacency_dictionary[vertex] = copy(structure.adjacency_dictionary[vertex])
            return possible_location
        # If the molecule does contain a subgraph which is isomorphic to the current CS
        # then check if this subgraph shares any vertices or edges with the structure
        for vertex in structure.adjacency_dictionary:
            neighbours_copy = copy(structure.adjacency_dictionary[vertex])
            for neighbour in structure.adjacency_dictionary[vertex]:
                if neighbour in self.path_structures[structure][molecule]:
                    if self.path_structures[structure][molecule][neighbour] in self.cs_locations[molecule]:
                        cs_neighbour = self.cs_locations[molecule][self.path_structures[structure][molecule][neighbour]]
                        neighbours_copy[cs_neighbour] = copy(structure.adjacency_dictionary[vertex][neighbour])
                        del neighbours_copy[neighbour]
            if self.path_structures[structure][molecule][vertex] in self.cs_locations[molecule]:
                cs_vertex = self.cs_locations[molecule][self.path_structures[structure][molecule][vertex]]
                possible_location.adjacency_dictionary[cs_vertex].update(neighbours_copy)
            elif self.path_structures[structure][molecule][vertex] not in self.cs_locations[molecule]:
                possible_location.vertex_to_graph(vertex)
                possible_location.adjacency_dictionary[vertex] = neighbours_copy
        return possible_location

    def _add_multiple_to_characteristic(self, structure):
        """
        Find locations to add a structure which appears multiple times in the molecules to characteristic substructure.

        :param structure: The structure which appears multiple times in the molecules and is to be added to CS
        :return: a list of graphs displaying the possible locations that the structure could have in the CS
        """
        print 'multiple'
        if len(self.multiple_vertices[structure]) < (len(self.molecules) * self.isomorphism_factor):
            return []
        repeats = set()
        for molecule in self.multiple_vertices[structure]:
            # Finds the number of different molecules associated  with one of the structure vertices
            # This indicates the number of instances of the structure that can be found in the molecule
            repeats.add(len(self.multiple_vertices[structure][molecule][structure.vertices()[0]]))
        repetitions = sorted(repeats, key=float)
        k_subgraphs = [structure]
        for k in repetitions:
            # A list of the molecules that contain exactly k incidences of the subgraph
            multi_molecules = [m for m in self.multiple_vertices[structure]
                               if len(self.multiple_vertices[structure][m].values()[0]) == k]
            for molecule in multi_molecules:
                # Create a graph representing the multiple copies of the structure
                path = (str(structure) + ' ') * k
                v = [vertex[index] for vertex in self.multiple_vertices[structure][molecule].values()
                     for index in range(len(self.multiple_vertices[structure][molecule].values()[0]))]
                vertices = []
                for vertex in v:
                    if vertex not in vertices:
                        vertices.append(vertex)
                multi_structure_tuple = self._create_structure(path, molecule, vertices)
                multi_structure = multi_structure_tuple[0]
                multi_vertices = multi_structure_tuple[1]
                unique_multi = self._check_structure_duplicates(multi_structure, molecule, multi_vertices)
                if unique_multi:
                    k_subgraphs.append(unique_multi)
        # Once there is a non repeated list of multiple structures made from structure try adding them to the CS
        possible_locations = {}
        for multi in k_subgraphs:
            for molecule in self.path_structures[multi]:
                # Create a possible location which is a combination of the CS and the multi_structure
                possible_location = self._add_single_to_characteristic(multi, molecule)
                possible_locations[possible_location] = multi
        if possible_locations:
            chosen_location = self._most_frequent_location(possible_locations.keys())
            self.characteristic_substructure = chosen_location
            self._add_cs_locations(possible_locations[chosen_location])

    def _add_cs_locations(self, structure):
        """
        Remembers the locations where molecules map to the characteristic substructure and the structure which is added.

        :param structure: the graph has been recently added to the characteristic substructure
        :return: None
        """
        self.cs_structures.append(structure)
        for molecule in self.path_structures[structure]:
            if molecule not in self.cs_locations:
                self.cs_locations[molecule] = {}
            for vertex in self.path_structures[structure][molecule]:
                # If this vertex of the structure has been added to the CS
                # then add the association to the cs_locations
                if vertex in self.characteristic_substructure.adjacency_dictionary.keys():
                    self.cs_locations[molecule][self.path_structures[structure][molecule][vertex]] = vertex

    def _most_frequent_location(self, possible_locations):
        """
        Finds the most frequent location of a structure given a list of graphs displaying possible locations.

        :param possible_locations: list of graphs which combine the CS and a structure
        :return: the graph which is chosen as most frequent
        """
        isomorphic_locations = {}
        for possible in possible_locations:
            for location in isomorphic_locations:
                if self._nx_isomorphism(possible, location):
                    isomorphic_locations[location] += 1
                    break
            else:
                isomorphic_locations[possible] = 1
        return OrderedDict(sorted(isomorphic_locations.items(), key=lambda x: x[1], reverse=True)).keys()[0]

    def _create_cs_results(self):
        """
        Creates the strings which will be used in the results output after calling find_characteristic_substructure.

        Uses the structures_results method to display the structures which are contained in the CS and the membership
        of all the molecules
        :return: list of display strings
        """
        string_list = ['Characteristic Substructure\n',
                       self._adjacency_dictionary_output(self.characteristic_substructure) + '\n',
                       'Structures which have been added to Characteristic Substructure\n']
        string_list.extend(self._structures_output(self.cs_structures))
        return string_list

    def _structures_output(self, structures):
        """
        Creates the strings which will be used in the results output after calling find_all_representative_structures

        Each molecule is presented with a bit array demonstrating the membership of the structures in the molecule
        :param structures: list of the representative structures
        :return: list of display strings
        """
        string_list = []
        counter = 0
        for structure in structures:
            string_list.append('Structure ' + str(counter) + '\n')
            string_list.append(self._adjacency_dictionary_output(structure) + '\n')
            counter += 1
        for molecule in self.molecules:
            string_list.append(str(molecule) + ': ')
            membership = []
            for structure in structures:
                if molecule in self.path_structures[structure]:
                    membership.append('1')
                else:
                    membership.append('0')
            string_list.append('(' + ''.join(membership) + ')' + '\n')
        return string_list

    def _adjacency_dictionary_output(self, structure):
        """
        Creates a string format of the adjacency dictionary of the given structure

        :param structure: molecule object which is the structure that is to be formatted
        :return: string representation of the molecule
        """
        string_list = []
        for vertex in structure.adjacency_dictionary:
            string_list.append(str(vertex) + ':\n')
            for neighbour in structure.neighbours(vertex):
                string_list.append('        ' + str(neighbour) + ': ' + str(structure.adjacency_dictionary[vertex][neighbour]) + '\n')
        display_string = ''.join(string_list)
        return display_string

    def _data_input(self):
        """
        Takes in a txt file, reads each line and stores the result in a list

        :return: a list of each of the lines of the text file
        """
        smiles_set = []
        reader = open(self.smiles_file, mode='rb')
        for line in reader:
            words = line.split(" ")
            smiles_set.append(words[0])
        reader.close()
        return smiles_set

    def _data_output(self, string_list, file_name):
        """
        Writes a list of strings to a file with the given filename

        :param string_list: list of display strings
        :return: None
        """
        display_string = ''.join(string_list)
        writer = open(file_name, mode='wb')
        writer.write(display_string)
        writer.close()


def argument_input():
    """
    Creates the algorithm object with the correct parameters and calls methods based on the command line arguments.

    The user can specify the smiles input file, whether they want to find characteristic substructure or all the
    representative structures (0 or 1) or they can specify the frequency threshold for paths and structures in molecules
    :return: None
    """
    c_structure = None
    all_structures = None
    if len(sys.argv) == 2:

        # Command: python algorithm.py 0/1
        if sys.argv[1] == 'C':
            cs = CSAlgorithm()
            c_structure = cs.find_characteristic_substructure()
        elif sys.argv[1] == 'R':
            cs = CSAlgorithm()
            all_structures = cs.find_all_representative_structures()

        # # Command: python algorithm.py 0.6
        # elif 0 < float(sys.argv[1]) < 1:
        #     cs = CSAlgorithm(threshold=float(sys.argv[1]))
        #     c_structure = cs.find_characteristic_substructure()
        #     all_structures = cs.find_all_representative_structures()

        # Command: python algorithm.py smiles_file
        else:
            cs = CSAlgorithm(smiles_file=sys.argv[1])
            c_structure = cs.find_characteristic_substructure()
            cs = CSAlgorithm(smiles_file=sys.argv[1])
            all_structures = cs.find_all_representative_structures()

    elif len(sys.argv) == 3:
        # Command: python algorithm.py 0/1 0.6
        # Want to use default input file but wish to choose threshold
        if sys.argv[1] == 'C' and 0 < float(sys.argv[2]) < 1:
            cs = CSAlgorithm(threshold=float(sys.argv[2]))
            c_structure = cs.find_characteristic_substructure()
        elif sys.argv[1] == 'R' and 0 < float(sys.argv[2]) < 1:
            cs = CSAlgorithm(threshold=float(sys.argv[2]))
            all_structures = cs.find_all_representative_structures()

        # Command: python algorithm.py smiles_file 0.6
        # Specify threshold but both characteristic substructure and representative structures are given
        elif sys.argv[1] != 'C' or sys.argv[1] != 'R' and 0 < float(sys.argv[2]) < 1:
            cs = CSAlgorithm(smiles_file=sys.argv[1], threshold=float(sys.argv[2]))
            c_structure = cs.find_characteristic_substructure()
            all_structures = cs.find_all_representative_structures()

        # Command: python algorithm.py smiles_file 0/1
        if sys.argv[2] == 'C':
            cs = CSAlgorithm(smiles_file=sys.argv[1])
            c_structure = cs.find_characteristic_substructure()
        elif sys.argv[2] == 'R':
            cs = CSAlgorithm(smiles_file=sys.argv[1])
            all_structures = cs.find_all_representative_structures()

    elif len(sys.argv) == 4:
        # Command: python algorithm.py smiles_file 0/1 0.6
        # O for characteristic substructure, 1 for representative structures
        cs = CSAlgorithm(smiles_file=sys.argv[1], threshold=float(sys.argv[3]))
        if sys.argv[2] == 'C':
            c_structure = cs.find_characteristic_substructure()
        elif sys.argv[2] == 'R':
            all_structures = cs.find_all_representative_structures()

    else:
        # Default smiles file, threshold and both commands are run
        cs = CSAlgorithm()
        c_structure = cs.find_characteristic_substructure()
        cs = CSAlgorithm()
        all_structures = cs.find_all_representative_structures()

    # Display the structures in the command line
    if c_structure:
        print 'Characteristic Substructure'
        print c_structure.adjacency_dictionary.keys()
    if all_structures:
        print 'Representative Structures'
        print all_structures

if __name__ == '__main__':
    argument_input()