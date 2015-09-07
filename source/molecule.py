import graph
from collections import deque


# Represents an atom within the molecule, includes the atoms chemical properties
class Atom(graph.Vertex):
    """
    A single atom of a molecule which includes the chemical information

    :param label: a string representing the atomic element
    :param isotope: the isotopic number if specified
    :param hydrogen: the number of explicit hydrogens
    :param charge: the charge of the atom
    :param aromatic: sets if the atom is part of an aromatic ring
    """
    def __init__(self, label, isotope=None, hydrogen=None, charge=None, aromatic=False):
        graph.Vertex.__init__(self, label)
        self.isotope = isotope
        self.hydrogen = hydrogen
        self.charge = charge
        self.aromatic = aromatic
        self.ring_break = False

    def __str__(self):
        """
        String representation based on atom type
        :return: string representation
        """
        if self.aromatic:
            return 'aromatic %s atom' % self.label
        else:
            return '%s atom' % self.label

    def __repr__(self):
        """
        Unique representation based on atom type
        :return: None
        """
        if self.aromatic:
            return '<Aromatic Atom %s at %s>' % (self.label, id(self))
        else:
            return '<Atom %s at %s>' % (self.label, id(self))


class Bond(graph.Edge):
    """
    A single bond of the molecule which has the bond type

    :param origin: one of the endpoints
    :param destination: one of the endpoints
    :param single: True if single bond
    :param double: True if double bond
    :param triple: True if triple bond
    :param quadruple: True if quadruple bond
    :param aromatic: True if aromatic bond
    """
    def __init__(self, origin, destination, single=False, double=False, triple=False, quadruple=False, aromatic=False):
        graph.Edge.__init__(self, origin, destination)
        self.single = single
        self.double = double
        self.triple = triple
        self.quadruple = quadruple
        self.aromatic = aromatic

    def __str__(self):
        """
        String representation that is different depending on the type of bond
        :return: string representation
        """
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
        """
        Unique representation that is different depending on the type of bond
        :return: unique representation of the bond
        """
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


class Molecule(graph.Graph):
    """
    A Molecule which inherits from the Graph class and contains atoms and bonds

    -smiles_string: the smiles string which is associated with the molecule
    """
    def __init__(self, smiles):
        graph.Graph.__init__(self)
        self._smiles_string = smiles

    # Return the original SMILES string
    def __str__(self):
        return self._smiles_string

    def __repr__(self):
        return '<Molecule %s at %s>' % (str(self), id(self))

    def add_atom(self, label, isotope=None, hydrogen=None, charge=None):
        """
        Adds a new atom to the molecule
        :param label: the atomic element
        :param isotope: the isotope number
        :param hydrogen: the number of explicit hydrogens
        :param charge: the charge of the atom
        :return: the new atom object that has been added
        """
        new_atom = Atom(label, isotope, hydrogen, charge)
        graph.Graph.vertex_to_graph(self, new_atom)
        return new_atom

    def add_aromatic_atom(self, label, isotope=None, hydrogen=None, charge=None):
        """
        Adds a new aromatic atom to the molecule
        :param label: the atomic element
        :param isotope: the isotope number
        :param hydrogen: the number of explicit hydrogens
        :param charge: the charge of the atom
        :return: the new aromatic atom
        """
        new_aromatic_atom = Atom(label, isotope, hydrogen, charge, aromatic=True)
        graph.Graph.vertex_to_graph(self, new_aromatic_atom)
        return new_aromatic_atom

    def add_single_bond(self, origin, destination):
        """
        Adds a single bond between the two provided atoms
        :param origin: an atom that is an endpoint
        :param destination: an atom that is an endpoint
        :return: the new bond object that has been added to the molecule
        """
        new_single = Bond(origin, destination, single=True)
        graph.Graph.edge_to_graph(self, origin, destination, new_single)
        return new_single

    def add_double_bond(self, origin, destination):
        """
        Adds a double bond between the two provided atoms
        :param origin: an atom that is an endpoint
        :param destination: an atom that is an endpoint
        :return: the new bond object that has been added to the molecule
        """
        new_double = Bond(origin, destination, double=True)
        graph.Graph.edge_to_graph(self, origin, destination, new_double)
        return new_double

    def add_triple_bond(self, origin, destination):
        """
        Adds a triple bond between the two provided atoms
        :param origin: an atom that is an endpoint
        :param destination: an atom that is an endpoint
        :return: the new bond object that has been added to the molecule
        """
        new_triple = Bond(origin, destination, triple=True)
        graph.Graph.edge_to_graph(self, origin, destination, new_triple)
        return new_triple

    def add_quadruple_bond(self, origin, destination):
        """
        Adds a quadruple bond between the two provided atoms
        :param origin: an atom that is an endpoint
        :param destination: an atom that is an endpoint
        :return: the new bond object that has been added to the molecule
        """
        new_quadruple = Bond(origin, destination, quadruple=True)
        graph.Graph.edge_to_graph(self, origin, destination, new_quadruple)
        return new_quadruple

    def add_aromatic_bond(self, origin, destination):
        """
        Adds an aromatic bond between the two provided atoms
        :param origin: an atom that is an endpoint
        :param destination: an atom that is an endpoint
        :return: the new bond object that has been added to the molecule
        """
        new_aromatic_bond = Bond(origin, destination, aromatic=True)
        graph.Graph.edge_to_graph(self, origin, destination, new_aromatic_bond)
        return new_aromatic_bond


