import graph
from collections import deque


# Represents an atom within the molecule, includes the atoms chemical properties
class Atom(graph.Vertex):
    def __init__(self, element, isotope=None, hydrogen=None, charge=None, aromatic=False):
        graph.Vertex.__init__(self, element)
        self.isotope = isotope
        self.hydrogen = hydrogen
        self.charge = charge
        self.ring_break = False
        self.aromatic = aromatic
    
    def __str__(self):
        if self.aromatic:
            return 'aromatic %s atom at position %s' % (self.element, self.position)
        else:
            return '%s atom at position %s' % (self.element, self.position)

    def __repr__(self):
        if self.aromatic:
            return '<Aromatic Atom %s at %s>' % (self.element, id(self))
        else:
            return '<Atom %s at %s>' % (self.element, id(self))


class Bond(graph.Edge):
    def __init__(self, origin, destination, single=False, double=False, triple=False, quadruple=False, aromatic=False):
        graph.Edge.__init__(self, origin, destination)
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


class Molecule(graph.Graph):
    def __init__(self, smiles):
        graph.Graph.__init__(self)
        self._smiles_string = smiles

    # Return the original SMILES string
    def __str__(self):
        return self._smiles_string

    def __repr__(self):
        return '<Molecule %s at %s>' % (str(self), id(self))

    def add_atom(self, element, isotope=None, hydrogen=None, charge=None):
        new_atom = Atom(element, isotope, hydrogen, charge)
        graph.Graph.vertex_to_graph(self, new_atom)
        return new_atom

    def add_aromatic_atom(self, element, isotope=None, hydrogen=None, charge=None):
        new_aromatic_atom = Atom(element, isotope, hydrogen, charge, aromatic=True)
        graph.Graph.vertex_to_graph(self, new_aromatic_atom)
        return new_aromatic_atom

    def add_single_bond(self, origin, destination):
        new_single = Bond(origin, destination, single=True)
        graph.Graph.edge_to_graph(self, origin, destination, new_single)
        return new_single

    def add_double_bond(self, origin, destination):
        new_double = Bond(origin, destination, double=True)
        graph.Graph.edge_to_graph(self, origin, destination, new_double)
        return new_double

    def add_triple_bond(self, origin, destination):
        new_triple = Bond(origin, destination, triple=True)
        graph.Graph.edge_to_graph(self, origin, destination, new_triple)
        return new_triple

    def add_quadruple_bond(self, origin, destination):
        new_quadruple = Bond(origin, destination, quadruple=True)
        graph.Graph.edge_to_graph(self, origin, destination, new_quadruple)
        return new_quadruple

    def add_aromatic_bond(self, origin, destination):
        new_aromatic_bond = Bond(origin, destination, aromatic=True)
        graph.Graph.edge_to_graph(self, origin, destination, new_aromatic_bond)
        return new_aromatic_bond

    def create_smiles(self):
        # Arbitrary vertex picked as the starting point
        start = self.adjacency_dictionary.keys()[0]
        parent = {}
        cycles = {}
        self._find_spanning_tree(start, parent, cycles)
        print parent
        print cycles
        # if len(parent) < self.size:
        #     parent = {}
        #     for u in self.adjacency_dictionary:
        #         if u not in parent:
        #             parent[u] = None
        #             self._find_spanning_tree(u, parent, cycles)
        # When vertex has a parent is None and isn't the start vertex then it is a dot break

    def _find_spanning_tree(self, s, parent, cycles):
        ring_number = 1
        parent[s] = None
        for v in self.vertices()[1:]:
            print 'starting here'
            print v
            if v not in parent:
                print 'deary'
                self._tree_visit(v, v, parent, cycles, ring_number)

    def _tree_visit(self, v, origin, parent, cycles, ring_number):
        for w in self.adjacency_dictionary[v]:
            print w
            print origin
            print v
            if w == origin and w != v:
                print 'cycle'
                if v in cycles:
                    string = cycles[v] + str(ring_number)
                    cycles[v] = string
                else:
                    cycles[v] = str(ring_number)
                if origin in cycles:
                    string = cycles[origin] + str(ring_number)
                    cycles[origin] = string
                else:
                    cycles[origin] = str(ring_number)
                ring_number += 1
            if w not in parent:
                print 'not a cycle'
                parent[w] = v
                self._tree_visit(w, origin, parent, cycles, ring_number)


    # # Based on Algorithms and Data Structures implementation of the DFS
    # def _depth_first_search(self, vertex, discovered, cycles, ring_number):
    #     for e in self.connecting_edges(vertex):
    #         v = e.opposite(vertex)
    #         if v not in discovered:
    # #             discovered[v] = e
    # #             self._depth_first_search(v, discovered, cycles, ring_number)
    # #         else:
    #             if v in cycles:
    #                 string = cycles[v] + str(ring_number)
    #                 cycles[v] = string
    #             else:
    #                 cycles[v] = str(ring_number)
    #             if e in cycles:
    #                 string = cycles[e] + str(ring_number)
    #                 cycles[e] = string
    #             else:
    #                 cycles[e] = str(ring_number)
    #             ring_number += 1


