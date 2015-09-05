from .graph import Vertex, Edge, Graph


# Represents an atom within the molecule, includes the atoms chemical properties
class Atom(Vertex):
    def __init__(self, element, isotope=None, hydrogen=None, charge=None, aromatic=False):
        Vertex.__init__(self, element)
        self.isotope = isotope
        self.hydrogen = hydrogen
        self.charge = charge
        self.ring_break = False
        self.aromatic = aromatic
    
    def __str__(self):
        if self.aromatic:
            return 'aromatic %s atom' % self.element
        else:
            return '%s atom' % self.element

    def __repr__(self):
        if self.aromatic:
            return '<Aromatic Atom %s at %s>' % (self.element, id(self))
        else:
            return '<Atom %s at %s>' % (self.element, id(self))


class Bond(Edge):
    def __init__(self, origin, destination, single=False, double=False, triple=False, quadruple=False, aromatic=False):
        Edge.__init__(self, origin, destination)
        self.single = single
        self.double = double
        self.triple = triple
        self.quadruple = quadruple
        self.aromatic = aromatic

    def __str__(self):
        if self.single:
            return 'single bond'
        elif self.double:
            return 'double bond'
        elif self.triple:
            return 'triple bond'
        elif self.quadruple:
            return 'quadruple bond'
        elif self.aromatic:
            return 'aromatic bond'

    def __repr__(self):
        if self.single:
            return '<Single Bond at %s>' % id(self)
        elif self.double:
            return '<Double Bond at %s>' % id(self)
        elif self.triple:
            return '<Triple Bond at %s>' % id(self)
        elif self.quadruple:
            return '<Quadruple Bond at %s>' % id(self)
        elif self.aromatic:
            return '<Aromatic Bond at %s>' % id(self)


class Molecule(Graph):
    def __init__(self, smiles):
        Graph.__init__(self)
        self._smiles_string = smiles

    # Return the original SMILES string
    def __str__(self):
        return self._smiles_string

    def __repr__(self):
        return '<Molecule %s at %s>' % (str(self), id(self))

    def add_atom(self, element, isotope=None, hydrogen=None, charge=None):
        new_atom = Atom(element, isotope, hydrogen, charge)
        Graph.vertex_to_graph(self, new_atom)
        return new_atom

    def add_aromatic_atom(self, element, isotope=None, hydrogen=None, charge=None):
        new_aromatic_atom = Atom(element, isotope, hydrogen, charge, aromatic=True)
        Graph.vertex_to_graph(self, new_aromatic_atom)
        return new_aromatic_atom

    def add_single_bond(self, origin, destination):
        new_single = Bond(origin, destination, single=True)
        Graph.edge_to_graph(self, origin, destination, new_single)
        return new_single

    def add_double_bond(self, origin, destination):
        new_double = Bond(origin, destination, double=True)
        Graph.edge_to_graph(self, origin, destination, new_double)
        return new_double

    def add_triple_bond(self, origin, destination):
        new_triple = Bond(origin, destination, triple=True)
        Graph.edge_to_graph(self, origin, destination, new_triple)
        return new_triple

    def add_quadruple_bond(self, origin, destination):
        new_quadruple = Bond(origin, destination, quadruple=True)
        Graph.edge_to_graph(self, origin, destination, new_quadruple)
        return new_quadruple

    def add_aromatic_bond(self, origin, destination):
        new_aromatic_bond = Bond(origin, destination, aromatic=True)
        Graph.edge_to_graph(self, origin, destination, new_aromatic_bond)
        return new_aromatic_bond


