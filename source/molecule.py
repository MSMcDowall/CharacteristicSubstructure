import graph


class Atom(graph.Vertex):
    def __init__(self, element, isotope=None, hydrogen=None, charge=None, ring_break=False):
        graph.Vertex.__init__(self, element)
        self._isotope = isotope
        self._hydrogen = hydrogen
        self._charge = charge
        self._ring_break = ring_break

    # Return the isotope number of the vertex
    @property
    def isotope(self):
        return self._isotope

    # Return the number of hydrogens attached to the vertex
    @property
    def hydrogen(self):
        return self._hydrogen

    # Return the charge of the element
    @property
    def charge(self):
        return self._charge

    # Return if atom was the break point of a ring
    @property
    def ring_break(self):
        return self.ring_break

    @ring_break.setter
    def ring_break(self, boolean):
        self._ring_break = boolean

    def __str__(self):
        return 'Atom with element %s' % self.element


class AromaticAtom(Atom):
    def __init__(self, element, isotope=None, hydrogen=None, charge=None):
        Atom.__init__(self, element, isotope, hydrogen, charge)

    def __str__(self):
        return 'Aromatic element with element %s' % self.element


class SingleBond(graph.Edge):
    def __init__(self, origin, destination):
        graph.Edge.__init__(self, origin, destination)

    def __str__(self):
        # return 'Single bond joining (%s) and (%s)' % (self._origin, self._destination)
        pass

class DoubleBond(graph.Edge):
    def __init__(self, origin, destination):
        graph.Edge.__init__(self, origin, destination)

        # Any methods particular to double bonds

class TripleBond(graph.Edge):
    def __init__(self, origin, destination):
        graph.Edge.__init__(origin, destination)

        # Any methods particular to triple bonds


class QuadrupleBond(graph.Edge):
    def __init__(self, origin, destination):
        graph.Edge.__init__(origin, destination)

        # Any methods particular to quadruple bonds


class AromaticBond(graph.Edge):
    def __init__(self, origin, destination):
        graph.Edge.__init__(self, origin, destination)

        # Any methods particular to aromatic bonds


class Molecule(graph.Graph):
    def __init__(self, string):
        graph.Graph.__init__(self)
        self._smiles_string = string        # The original SMILES string

    # Return the original SMILES string
    @property
    def smiles_string(self):
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
        new_aromatic_bond = AromaticAtom(origin, destination)
        graph.Graph.edge_to_graph(self, origin, destination, new_aromatic_bond)
        return new_aromatic_bond

