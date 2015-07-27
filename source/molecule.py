import graph


# Represents an atom within the molecule, includes the atoms chemical properties
class Atom(graph.Vertex):
    def __init__(self, element, isotope=None, hydrogen=None, charge=None):
        graph.Vertex.__init__(self, element)
        self.isotope = isotope
        self.hydrogen = hydrogen
        self.charge = charge
        self.ring_break = False

    # Return if atom was the break point of a ring
    @property
    def ring_break(self):
        return self.ring_break

    @ring_break.setter
    def ring_break(self, boolean):
        self._ring_break = boolean
    
    def __str__(self):
        return 'Atom element %s at position %s' % (self.element, self.position)


# Represents an aromatic atom and inherits the chemical properties that are present in an atom
class AromaticAtom(Atom):
    def __init__(self, element, isotope=None, hydrogen=None, charge=None):
        Atom.__init__(self, element, isotope, hydrogen, charge)

    def __str__(self):
        return 'Aromatic element %s at position %s' % (self.element, self.position)


class SingleBond(graph.Edge):
    def __init__(self, origin, destination):
        graph.Edge.__init__(self, origin, destination)

    def __str__(self):
        return 'Single bond ' + graph.Edge.__str__(self)

class DoubleBond(graph.Edge):
    def __init__(self, origin, destination):
        graph.Edge.__init__(self, origin, destination)

    def __str__(self):
        return 'Double bond ' + graph.Edge.__str__(self)

class TripleBond(graph.Edge):
    def __init__(self, origin, destination):
        graph.Edge.__init__(self, origin, destination)

    def __str__(self):
        return 'Triple bond ' + graph.Edge.__str__(self)


class QuadrupleBond(graph.Edge):
    def __init__(self, origin, destination):
        graph.Edge.__init__(self, origin, destination)

    def __str__(self):
        return 'Quadruple bond ' + graph.Edge.__str__(self)


class AromaticBond(graph.Edge):
    def __init__(self, origin, destination):
        graph.Edge.__init__(self, origin, destination)

    def __str__(self):
        return 'Aromatic bond ' + graph.Edge.__str__(self)


class Molecule(graph.Graph):
    def __init__(self, smiles):
        graph.Graph.__init__(self)
        self._smiles_string = smiles        # The original SMILES string

    # Return the original SMILES string
    def __str__(self):
        return self._smiles_string

    def add_atom(self, element, isotope=None, hydrogen=None, charge=None):
        new_atom = Atom(element, isotope, hydrogen, charge)
        graph.Graph.vertex_to_graph(self, new_atom)
        return new_atom

    def add_aromatic_atom(self, element, isotope=None, hydrogen=None, charge=None):
        new_aromatic_atom = AromaticAtom(element, isotope, hydrogen, charge)
        graph.Graph.vertex_to_graph(self, new_aromatic_atom)
        return new_aromatic_atom

    def add_single_bond(self, origin, destination):
        new_single = SingleBond(origin, destination)
        graph.Graph.edge_to_graph(self, origin, destination, new_single)
        return new_single

    def add_double_bond(self, origin, destination):
        new_double = DoubleBond(origin, destination)
        graph.Graph.edge_to_graph(self, origin, destination, new_double)
        return new_double

    def add_triple_bond(self, origin, destination):
        new_triple = TripleBond(origin, destination)
        graph.Graph.edge_to_graph(self, origin, destination, new_triple)
        return new_triple

    def add_quadruple_bond(self, origin, destination):
        new_quadruple = QuadrupleBond(origin, destination)
        graph.Graph.edge_to_graph(self, origin, destination, new_quadruple)
        return new_quadruple

    def add_aromatic_bond(self, origin, destination):
        new_aromatic_bond = AromaticBond(origin, destination)
        graph.Graph.edge_to_graph(self, origin, destination, new_aromatic_bond)
        return new_aromatic_bond

