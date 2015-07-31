# coding=utf-8
from smiles_parser import Parser
import molecule
import subprocess
from collections import OrderedDict
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
        self.molecules = []                 # All the given molecules
        self.paths = {}                     # All the paths from all the molecules with their lengths
        self.path_structures = {}           # Dictionary of path structures which has record of isomorphic molecule
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
        structure = molecule.Molecule(path)
        old_molecule_map = {}           # Uses the original vertices as keys, structure vertices as values
        new_molecule_map = {}           # Uses the structure vertices as keys, original vertices as values
        for atom in vertices:
            if isinstance(atom, molecule.Atom):
                new_atom = structure.add_atom(atom.element)
                old_molecule_map[atom] = new_atom
                new_molecule_map[new_atom] = atom
            elif isinstance(atom, molecule.AromaticAtom):
                new_atom = structure.add_aromatic_atom(atom.element)
                old_molecule_map[atom] = new_atom
                new_molecule_map[new_atom] = atom
        for atom in vertices:
            for neighbour in vertices:
                edge = mole.contains_edge(atom, neighbour)      # Test if there is an edge to other atoms in path
                if isinstance(edge, molecule.SingleBond):
                    structure.add_single_bond(old_molecule_map[atom], old_molecule_map[neighbour])
                elif isinstance(edge, molecule.DoubleBond):
                    structure.add_double_bond(old_molecule_map[atom], old_molecule_map[neighbour])
                elif isinstance(edge, molecule.TripleBond):
                    structure.add_triple_bond(old_molecule_map[atom], old_molecule_map[neighbour])
                elif isinstance(edge, molecule.QuadrupleBond):
                    structure.add_quadruple_bond(old_molecule_map[atom], old_molecule_map[neighbour])
                elif isinstance(edge, molecule.AromaticBond):
                    structure.add_aromatic_bond(old_molecule_map[atom], old_molecule_map[neighbour])
        return structure, new_molecule_map

    def create_multiple_structures(self):
        # multiple_structures = {structure: mole: [[strut, vert map]...]}
        for structure in self.multiple_structures:
            if len(self.multiple_structures[structure]) > len(self.molecules)*self.isomorphism_factor:
                pass

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
                # Does the isomorphic structure already contain a reference to the current molecule?
                # This would mean that the molecule has multiple occurrences of the same structure
                # if mole in self.path_structures[structure]:
                #     print 'this structure in molecule 2 times'
                #     if structure in self.multiple_structures and mole in self.multiple_structures[structure]:
                #         self.multiple_structures[structure][mole].append([pattern, vertices])
                #     else:
                #         self.multiple_structures[structure] = {mole: [[pattern, vertices]]}
                #     print type(self.multiple_structures[structure])
                #     print type(self.multiple_structures[structure][mole])
                # else:
                #     print 'this structure not before in molecule'
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
                print 'NEW MOLECULE'
                path_vertices = [pair[1] for pair in mole.paths if pair[0] == path]    # Lists of path vertices
                print path_vertices
                for vertices in path_vertices:
                    structure_tuple = self.create_structure(path, mole, vertices)
                    structure = structure_tuple[0]
                    vertices = structure_tuple[1]
                    # Test if the structure has already been encountered - change here between LAD and nx isomorphism
                    self.create_nx_graph(structure)
                    self.nx_isomorphism(structure, mole, vertices)
        if self.multiple_structures:
            print 'multiple thingys'
            print self.multiple_structures
        for structure in self.path_structures:
            relative_frequency = len(self.path_structures[structure].keys())/float(len(self.molecules))
            if float(relative_frequency) >= self.threshold:
                rep_structures[structure] = relative_frequency
        # Store relative frequency of each structure (induced by path in rep_paths) as value in dictionary
        # Sort dictionary based on frequency highest to lowest
        representative_structures = OrderedDict(sorted(rep_structures.items(), key=lambda x: x[1], reverse=True))
        return representative_structures

    def add_structure_to_characteristic(self, sorted_list):
        # Use structures stored in order of frequency to decide where it should be added to characteristic substructure
        # When a structure is added to the CS, any molecule that is sub isomorphic to it will have its vertices swapped
        # These molecules will instead contain CS vertices
        for structure in sorted_list:
            for mole in self.path_structures[structure]:
                # print self.path_structures[structure][mole]
                # vertex_map = self.path_structures[structure][mole]
                # for key in vertex_map:
                #     mole.swap_vertex(vertex_map[key], key)
                pass

    def find_characteristic_substructure(self):
        smiles_set = []
        reader = open('SMILES', mode='rb')
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
                self.add_structure_to_characteristic(sorted_list)
                length -= self.step
            else:
                length -= 1

if __name__ == '__main__':
    path_finder = CharacteristicSubstructure(threshold=0.3, length_start=3, length_end=3, step=1)
    path_finder.find_characteristic_substructure()


