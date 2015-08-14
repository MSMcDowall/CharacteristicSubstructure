# coding=utf-8
from smiles_parser import Parser
import molecule as m
from collections import OrderedDict, Counter
from copy import copy, deepcopy
import networkx as nx
import networkx.algorithms.isomorphism as iso
#from draw_molecule import draw_molecule as draw


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
        # All the given molecules
        self.molecules = [] 
        # All the paths from all the molecules with their lengths                
        self.paths = {} 
        # Dictionary of path structures which has record of isomorphic molecules and mapping from structure vertices to molecule vertices 
        # {path structure: {molecule: {structure vertex: molecule vertex}}}
        self.path_structures = {}  
        # Dictionary of the structures which appear multiple times in molecules and a list of the molecule vertices that map to the vertices 
        # {structure: {molecule: {structure vertex : [molecule vertices]}}}        
        self.multiple_structures = {}       
        # NetworkX graph version of the path structure, path structure is key
        self.structure_nx = {}              

    def find_graphs_paths(self, smiles_set):
        """
        For each SMILES string a molecule object is created and all of its paths found
        
        :param smiles_set: a list of SMILES strings
        :return: a dictionary containing all of the paths as strings with their length as value
        """
        for smiles in smiles_set:
            mole = Parser().parse_smiles(smiles)
            #draw(mole)
            self.molecules.append(mole)
            path_dict = mole.find_all_paths()
            self.paths.update(path_dict)
        return self.paths

    def find_representative_paths(self, length):
        """
        Find the paths which occur with a high enough frequency
        
        :param length: an integer which sets the length of the paths that should be considered
        :return: a list of the paths as strings which are representative 
        """
        representative_paths = []
        for path in self.paths:
            # Check all the paths which are of the chosen length
            if self.paths[path] == length:      
                counter = 0
                for mole in self.molecules:
                    # Search the tuples for each path which consists of path string and vertices present in path
                    if path in (pair[0] for pair in mole.paths):        
                        counter += 1
                if float(counter)/len(self.molecules) >= self.threshold:
                    representative_paths.append(path)
        return representative_paths

    def find_representative_structures(self, rep_paths):
        """
        Find the path structures which appear with a high enough frequency 
        Ensures there are no duplicates in the dictionary of path structures so that the frequency count is accurate
        
        :param rep_paths: a list of the paths which are representive
        :return: a list of path structures which are representative sorted into order of frequency, highest to lowest
        """
        print 'rep_paths'
        print rep_paths
        rep_structures = {}
        for path in rep_paths:
            for mole in self.molecules:
                print 'mole to be structured'
                print str(mole)
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
            print self.multiple_structures
            # self._create_multiple_structures
        for structure in self.path_structures:
            relative_frequency = len(self.path_structures[structure].keys())/float(len(self.molecules))
            if float(relative_frequency) >= self.threshold:
                # Store relative frequency of each structure (induced by path in rep_paths) as value in dictionary
                rep_structures[structure] = relative_frequency
        # Sort dictionary based on frequency highest to lowest and return the path structures only
        return OrderedDict(sorted(rep_structures.items(), key=lambda x: x[1], reverse=True)).keys()

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
        self.find_graphs_paths(smiles_set)
        length = self.length_start
        while length >= self.length_end:
            representative_paths = self.find_representative_paths(length)
            print 'found representative paths'
            sorted_list = self.find_representative_structures(representative_paths)
            # After considering paths of this length test to see if there are representative substructures
            # If there are no rep structures then decrease stepwise, if there is increase the step size
            if sorted_list:
                print 'sorted list of representative structures'
                self.characteristic_substructure = copy(sorted_list[0])
                print 'lets swap CS'
                self._swap_molecule_vertices(sorted_list[0])
                for structure in sorted_list[1:]:
                    print 'adding more to CS'
                    self._add_structure_to_characteristic(structure)
                    print 'edited new cs struct'
                    self._swap_molecule_vertices(structure)
                length -= self.step
            else:
                length -= 1
        return self.characteristic_substructure
        
    
    # TODO restructure to check for my multiples and then go to either of 2 seperate methods
    # Leave swapping and location comparison in place 
    def _add_structure_to_characteristic(self, structure):
        """
        Adds the given structure to the characteristic substructure in the location where it appears most frequently in the molecules
        A different method is called depending on if the structure appears once or many times in the molecules
        A list of the possible locations on the CS where the structure could be added is created
        This list is then sorted in terms of frequency and the CS with the most common position of the structure becomes the new CS
        
        :param structure: the molecule object which is to be added to the characteristic substructure
        """
        possible_locations = []    # List of possible locations of subgraphs
        for mole in self.path_structures[structure]:
            # Create copy of the structure so that changes can be made without altering original structure
            structure_copy = m.Molecule(str(structure))
            for vertex in structure.adjacency_dictionary:
                structure_copy.adjacency_dictionary[vertex] = {}
                for neighbour in structure.adjacency_dictionary[vertex]:
                    structure_copy.adjacency_dictionary[vertex][neighbour] = copy(structure.adjacency_dictionary[vertex][neighbour])
            # Does structure vertices map to any CS vertices?
            # If they do then change them to the CS vertices
            print 'creating copy of structure'
            for vertex in self.path_structures[structure][mole]:
                if self.path_structures[structure][mole][vertex] in self.characteristic_substructure.adjacency_dictionary:
                    # SWAP VERTEX ERROR HERE for multiple structure
                    structure_copy.swap_vertex(vertex, self.path_structures[structure][mole][vertex])
            # Create a copy of the characteristic substructure which will have the structure added to it
            print 'creating location'
            possible_location = m.Molecule(str(self.characteristic_substructure))
            for vertex in self.characteristic_substructure.adjacency_dictionary:
                possible_location.adjacency_dictionary[vertex] = {}
                for neighbour in self.characteristic_substructure.adjacency_dictionary[vertex]:
                    possible_location.adjacency_dictionary[vertex][neighbour] = copy(self.characteristic_substructure.adjacency_dictionary[vertex][neighbour])
            # For vertices in structure not in CS
            # Add to CS also append bonds for structure vertices
            for vertex in structure_copy.adjacency_dictionary:
                if vertex in self.characteristic_substructure.adjacency_dictionary:
                    possible_location.adjacency_dictionary[vertex].update(structure_copy.adjacency_dictionary[vertex])
                else:
                    possible_location.adjacency_dictionary[vertex] = structure_copy.adjacency_dictionary[vertex]
            possible_locations.append(possible_location)
        for possible in possible_locations:
            self._create_nx_graph(possible)
        isomorphic_locations = {}
        for possible in possible_locations:
            for location in isomorphic_locations:
                if self._nx_isomorphism(possible, location):
                    isomorphic_locations[location] += 1
                    break
            else:
                isomorphic_locations[possible] = 1
        sorted_locations = OrderedDict(sorted(isomorphic_locations.items(), key=lambda x: x[1], reverse=True)).keys()
        self.characteristic_substructure = sorted_locations[0]


    def _create_structure(self, path, mole, vertices):
        """
        Creates a path strucuture which is a molecule object representing the form that the given path takes within a molecule
        
        :param path: the string which represents the path (Entered as the SMILES string for the path structure)
        :param molecule: the molecule object which contains the given path and vertices
        :param vertices: a list of vertices which are encountered in the path of the molecule
        :return: a tuple of the newly created path structure and a dictionary mapping the structures vertices to the molecule vertices 
        """
        print 'into create structure'
        structure = m.Molecule(path)
        # Uses the original vertices as keys, structure vertices as values, used to create the edges
        old_molecule_map = {}  
        # Uses the structure vertices as keys, original vertices as values, used to find the molecule vertices if the structure is given         
        new_molecule_map = {}           
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
                # Test if there is an edge to any of the other atoms in the path
                edge = mole.contains_edge(atom, neighbour)
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

    # TODO Put in place while adding to CS rather than before
    def _create_multiple_structures(self):
        # multiple_structures = {structure: mole: [[strut, vert map]...]}
        for structure in self.multiple_structures:
            # "If at least k>= 1 isomorphic subgraphs within the same molecule graph in at least |M|.iso graphs occur:"
            # If the subgraph structure appears in |M|.iso graphs then multiple structures will be created
            if len(self.path_structures[structure]) > len(self.molecules)*self.isomorphism_factor:
                for mole in self.multiple_structures[structure]:
                    # TODO Get longest set as value of vertex
                    k = len(self.multiple_structures[structure][mole])
                    path = (str(structure) + ' ') * k
                    new_structure = m.Molecule(path)
                    for vertex in structure.adjacency_dictionary:
                        new_structure.adjacency_dictionary[vertex] = {}
                        for neighbour in structure.adjacency_dictionary[vertex]:
                            new_structure.adjacency_dictionary[vertex][neighbour] = copy(structure.adjacency_dictionary[vertex][neighbour])
                    new_structure.size += structure.size
                    vertices_mapping = self.path_structures[structure][mole].copy()
                    # TODO Change to take in to account the new multiple structures setup
                    for pair in self.multiple_structures[structure][mole]:
                        substructure = pair[0]
                        # Ensure positions in the structure are unique so they can be used with the nx graph
                        for vertex in substructure.adjacency_dictionary:
                            vertex.position += new_structure.size
                            new_structure.adjacency_dictionary[vertex] = {}
                            for neighbour in substructure.adjacency_dictionary[vertex]:
                                new_structure.adjacency_dictionary[vertex][neighbour] = copy(substructure.adjacency_dictionary[vertex][neighbour])
                        new_structure.size += substructure.size
                        vertices_mapping.update(pair[1])
                    draw(new_structure)
                    # Remove that mole entry from path structures dictionary so it is not repeated
                    # del self.path_structures[structure][mole]
                    # Add the new larger structure into the dictionary
                    self._create_nx_graph(new_structure)
                    # self._check_structure_duplicates(new_structure, mole, vertices_mapping)
                    # draw(new_structure)

    def _create_nx_graph(self, path_structure):
        """
        Creates a copy of the molecule object as an NetworkX graph
        The position attirbute of each structure vertices is used as the index of the NetworkX graph
        
        :param path_structure: the molecule object which will be turned into a NetworkX graph
        """
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

    
    def _nx_isomorphism(self, pattern, target):
        """
        Uses the NetworkX isomorphism algorithm to check if the pattern graph and the target graph are isomorphic
        The faster_could_be_isomorphic method is used to discount two structures if they could not possibly be isomorphic
        
        :param pattern: a molecule object which is to be tested for isomorphism
        :param target: a molecule object which pattern graph is to be compared against
        :return: None if the graphs are not isomorphic
        :return: a dictionary which maps the indices of the two NetworkX graphs together if they are isomorphic
        """
        # Error calling without making structure
        if not nx.faster_could_be_isomorphic(self.structure_nx[pattern], self.structure_nx[target]):
                # Graphs are definitely not isomorphic
                return None
        matcher = iso.GraphMatcher(self.structure_nx[pattern], self.structure_nx[target],
                                   node_match=iso.categorical_node_match('element', 'C'),
                                   edge_match=iso.categorical_edge_match('type', 'single'))
        if matcher.is_isomorphic():
            return matcher.mapping

    def _check_structure_duplicates(self, pattern, mole, vertices):
        temporary_structure_dict = self.path_structures.copy()
        print 'structures'
        print self.path_structures
        for structure in self.path_structures.keys():
            print 'back again'
            isomorphic_mapping = self._nx_isomorphism(pattern, structure)
            if not isomorphic_mapping:
                # Continues to next structure for testing as they are not isomorphic
                continue
            if isomorphic_mapping:
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
                print 'altered mapping'
                print altered_mapping
                if mole in temporary_structure_dict[structure]:
                    # If the molecule vertices are different to the vertices currently in dictionary
                    # then there are multiple subgraphs in molecule isomorphic to pattern
                    if Counter(vertices.values()) != Counter(temporary_structure_dict[structure][mole].values()):
                        print 'matching ness'
                        print pattern.adjacency_dictionary
                        print structure.adjacency_dictionary
                        if structure in self.multiple_structures and mole in self.multiple_structures[structure]:
                            print 'multiple srtucture'
                            print self.multiple_structures[structure][mole]
                            for key in altered_mapping:
                                print self.multiple_structures[structure][mole][key]
                                self.multiple_structures[structure][mole][key].add(altered_mapping[key])
                        else:
                            print 'new multiple structure'
                            # One extra copy of the structure is being added so TEMPORARY solution: make them sets
                            if structure in self.multiple_structures:
                                self.multiple_structures[structure][mole] = {}
                            elif structure not in self.multiple_structures:
                                self.multiple_structures[structure] = {mole: {}}
                            for key in altered_mapping:
                                self.multiple_structures[structure][mole][key] = set([altered_mapping[key],
                                                                self.path_structures[structure][mole][key]])
                            print 'new structy ness'
                            print self.multiple_structures
                else:
                    temporary_structure_dict[structure][mole] = altered_mapping
                break
        else:
            temporary_structure_dict[pattern] = {mole: vertices}
        self.path_structures = temporary_structure_dict.copy()

    def _data_input(self):
        """
        Takes in a txt file, reads each line and stores the result in a list
        
        :return: a list of each of the lines of the text file
        """
        smiles_set = []
        reader = open('SMILES.txt', mode='rb')
        for line in reader:
            smiles_set.append(line.rstrip())
        reader.close()
        return smiles_set

    def _data_output(self):
        pass

    def _swap_molecule_vertices(self, structure):
        """
        Changes the vertices of the molecules with the vertices of the characteristic substructure that they map to
        
        :param structure: the graph object that has been added to the CS
        """
        for mole in self.path_structures[structure]:
            for key in self.path_structures[structure][mole]:
                # If this vertex of the structure has been added to the CS
                # then change the molecule that are associated with it
                if key in self.characteristic_substructure.adjacency_dictionary.keys():
                    # Change the vertices in the molecule
                    mole.swap_vertex(self.path_structures[structure][mole][key], key)
                    self._swap_path_structure(self.path_structures[structure][mole][key], key)


    def _swap_path_structure(self, old, new):
        """
        Changes each occasion of the molecule vertex with the CS vertex that has replaced it
        """
        for structure in self.path_structures:
            for molecule in self.path_structures[structure]:
                for vertex in self.path_structures[structure][molecule]:
                    if self.path_structures[structure][molecule][vertex] == old:
                        self.path_structures[structure][molecule][vertex] = new
        for structure in self.multiple_structures:
            for molecule in self.multiple_structures[structure]:
                for vertex in self.multiple_structures[structure][molecule]:
                    if old in self.multiple_structures[structure][molecule][vertex]:
                        self.multiple_structures[structure][molecule][vertex].add(new)
                        self.multiple_structures[structure][molecule][vertex].remove(old)

if __name__ == '__main__':
    path_finder = CharacteristicSubstructure(threshold=0.3, length_start=3, length_end=3)
    struct = path_finder.find_characteristic_substructure()
    print struct
