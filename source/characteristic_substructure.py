# coding=utf-8
from parser import Parser
import molecule
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
        self.path_structures = {}           # Dictionary of path structures which has a dictionary of molecules that it
                                            # appears in as its value

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
        print 'molecule map'
        print molecule_map
        for atom in vertices:
            for neighbour in vertices:
                edge = mole.contains_edge(atom, neighbour)
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
        print structure
        for vertex in structure.vertices:
            print vertex
        self.path_structures[structure] = {mole: vertices}
        return structure

    def multiple_subgraphs(self, path, mole, positions):
        pass

    def find_representative_structures(self, rep_paths):
        resresentative_structures = {}
        for path in rep_paths:
            for mole in self.molecules:
                #if path in (pair[0] for pair in mole.paths):
                    path_vertices = [pair[1] for pair in mole.paths if pair[0] == path]    # Lists of path vertices
                    # If path_vertices > 1 go to special method
                    if len(path_vertices) > 1:
                        new_structure = self.multiple_subgraphs(path, mole, path_vertices)
                    elif len(path_vertices) > 0:
                        new_structure = self.create_structure(path, mole, path_vertices[0])
                    # If new_structure is isomorphic to another subgraph in path_structures.keys()
                    #
        # Create structures from path
        # Test for subgraph in each molecule
        # Store location of each substructure that is isomorphic to representative structure as:
        # Dictionary within dictionary {path_structure: {molecule : position}}
        # Consider here what happens if there is more than one subgraph isomorphic to path structure in molecule
        # Store relative frequency of each structure (induced by path in rep_paths) as value in dictionary 'rep_strut'
        # Sort 'rep_strut' based on frequency highest to lowest
        # OrderedDict(sorted(rep_strut.items(), key=lambda x: x[1], reverse=True))
        return resresentative_structures     # Return the sorted collection of path structures

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
    path_finder = CharacteristicSubstructure(path_threshold=0.5)
    print path_finder.find_graphs_paths(['CNO', 'CPOS', 'ClFIN', 'POSN'])
    for mole in path_finder.molecules:
        print mole
        draw(mole)
    rep_paths = path_finder.find_representative_paths(3)
    print rep_paths
    structures = path_finder.find_representative_structures(rep_paths)
    for subgraph in path_finder.path_structures.keys():
        draw(subgraph)
