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
        self.length_start = length_start
        self.length_end = length_end
        self.step = step
        self.threshold = threshold
        self.isomorphism_factor = isomorphism_factor
        self.characteristic_substructure = None
        self.molecules = []                 # All the given molecules
        self.paths = {}                     # All the paths from all the molecules with their lengths
        self.path_structures = {}           # Dictionary of path structures which has record of isomorphic molecules
        self.isomorphic_structures = {}     # Dictionary of all the structures which are isomorphic to one another
        self.structure_text = {}            # Text version of graph for LAD isomorphism
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

    def create_structure(self, path, mole, vertices):
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

    def create_multiple_structures(self):
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
                        new_structure.adjacency_dictionary.update(substructure.adjacency_dictionary)
                        new_structure.size += substructure.size
                        vertices_mapping.update(pair[1])
                    # Remove that mole entry from path structures dictionary so it is not repeated
                    del self.path_structures[structure][mole]
                    # Add the new larger structure into the dictionary
                    self.create_nx_graph(new_structure)
                    self.check_structure_duplicates(new_structure, mole, vertices_mapping)
                    # draw(new_structure)

    def create_nx_graph(self, path_structure):
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

    def nx_isomorphism(self,pattern, target):
        if not nx.faster_could_be_isomorphic(pattern, target):
                # Graphs are definitely not isomorphic
                return False
        matcher = iso.GraphMatcher(pattern, target,
                                   node_match=iso.categorical_node_match('element', 'C'),
                                   edge_match=iso.categorical_edge_match('type', 'single'))
        if matcher.is_isomorphic():
            return True
        else:
            return False

    def check_structure_duplicates(self, pattern, mole, vertices):
        temporary_structure_dict = self.path_structures.copy()
        for structure in self.path_structures.keys():
            isomorphic = self.nx_isomorphism(self.structure_nx[pattern], self.structure_nx[structure])
            if not isomorphic:
                continue    # Continues to next structure for testing as they are not isomorphic
            if isomorphic:
                if mole in temporary_structure_dict[structure]:
                    if Counter(vertices.values()) == Counter(temporary_structure_dict[structure][mole].values()):
                        print 'same vertices'
                        temporary_structure_dict[structure][mole] = vertices
                    else:
                        print 'different vertices'
                        if structure in self.multiple_structures and mole in self.multiple_structures[structure]:
                            self.multiple_structures[structure][mole].append([pattern, vertices])
                        else:
                            self.multiple_structures[structure] = {mole: [[pattern, vertices]]}
                else:
                    temporary_structure_dict[structure][mole] = vertices
                break
        else:
            temporary_structure_dict[pattern] = {mole: vertices}
        self.path_structures = temporary_structure_dict.copy()

    def find_representative_structures(self, rep_paths):
        rep_structures = {}
        for path in rep_paths:
            print path
            for mole in self.molecules:
                path_vertices = [pair[1] for pair in mole.paths if pair[0] == path]    # Lists of path vertices
                for vertices in path_vertices:
                    structure_tuple = self.create_structure(path, mole, vertices)
                    structure = structure_tuple[0]
                    print 'new structure'
                    print mole
                    print structure
                    print structure.adjacency_dictionary.keys()
                    vertices = structure_tuple[1]
                    # Test if the structure has already been encountered - change here between LAD and nx isomorphism
                    self.create_nx_graph(structure)
                    self.check_structure_duplicates(structure, mole, vertices)
        if self.multiple_structures:
            self.create_multiple_structures()
        for structure in self.path_structures:
            print 'check numbers'
            print structure
            print len(self.path_structures[structure].keys())
            relative_frequency = len(self.path_structures[structure].keys())/float(len(self.molecules))
            if float(relative_frequency) >= self.threshold:
                rep_structures[structure] = relative_frequency
        # Store relative frequency of each structure (induced by path in rep_paths) as value in dictionary
        # Sort dictionary based on frequency highest to lowest
        representative_structures = OrderedDict(sorted(rep_structures.items(), key=lambda x: x[1], reverse=True))
        print representative_structures
        return representative_structures

    def swap_path_structure(self, old, new):
        for structure in self.path_structures:
            for molecule in self.path_structures[structure]:
                for vertex in self.path_structures[structure][molecule]:
                    if self.path_structures[structure][molecule][vertex] == old:
                        self.path_structures[structure][molecule][vertex] = new

    def add_structure_to_characteristic(self, structure):
        # When a structure is added to the CS, any molecule that is sub isomorphic to it will have its vertices swapped
        # These molecules will instead contain CS vertices
        print 'strcuty struct'
        print self.path_structures[structure]
        print structure.adjacency_dictionary.keys()
        if not self.characteristic_substructure:
            self.characteristic_substructure = copy(structure)
            print 'CS made'
            print self.characteristic_substructure.adjacency_dictionary
        # elif self.characteristic_substructure:
        #     print 'CHARACTERISTIC SUBSTRUCTURE'
        #     possible_locations = {}    # List of possible locations of subgraphs
        #     for mole in self.path_structures[structure]:
        #         cs_vertices = [vertex for vertex in mole.adjacency_dictionary.keys()    # Overlap between molecule & CS
        #                        if vertex in self.characteristic_substructure.adjacency_dictionary.keys()]
        #         possible_location = copy(self.characteristic_substructure)
        #         print possible_location.adjacency_dictionary
                # Does structure vertices contain any from cs vertices?
                # For vertices in structure not in CS
                # Add to CS also append bonds for structure vertices
        # Compare the frequency of possible locations
        # Isomorphism algorithm
        print 'CS vertices'
        print self.characteristic_substructure.adjacency_dictionary.keys()
        for mole in self.path_structures[structure]:
            print 'path struct'
            print self.path_structures[structure][mole]
            for key in self.path_structures[structure][mole]:
                print repr(key)
                # If this vertex of the structure has been added to the CS
                if key in self.characteristic_substructure.adjacency_dictionary.keys():
                    print 'changing vertices'
                    print self.path_structures[structure][mole][key]
                    # Change the vertices in the molecule
                    mole.swap_vertex(self.path_structures[structure][mole][key], key)
                    self.swap_path_structure(self.path_structures[structure][mole][key], key)
                    print self.path_structures[structure][mole][key]

    def find_characteristic_substructure(self):
        smiles_set = []
        reader = open('SMILES.txt', mode='rb')
        for line in reader:
            smiles_set.append(line.rstrip())
        reader.close()
        self.find_graphs_paths(smiles_set)
        length = self.length_start
        while length >= self.length_end:
            representative_paths = self.find_representative_paths(length)
            sorted_list = self.find_representative_structures(representative_paths)
            # After considering paths of this length test to see if there are representative substructures
            # If there are no rep structures then decrease stepwise, if there is increase the step size
            if sorted_list:
                print 'list'
                print sorted_list
                for structure in sorted_list:
                    print 'sorted'
                    print structure
                    self.add_structure_to_characteristic(structure)
                length -= self.step
            else:
                length -= 1

if __name__ == '__main__':
    path_finder = CharacteristicSubstructure(threshold=0.3, length_start=3, length_end=3, step=1)
    path_finder.find_characteristic_substructure()


