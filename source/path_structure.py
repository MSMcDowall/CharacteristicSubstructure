import molecule as m
import smiles_parser
import networkx as nx


class PathStructure(object):
    def __init__(self, structure):
        self.structure = structure
        self.molecules = {}
        self.nx_graph = None

    def create_nx_graph(self):
        """
        Creates a copy of the molecule object as an NetworkX graph
        The position attirbute of each structure vertices is used as the index of the NetworkX graph

        :param structure: the molecule object which will be turned into a NetworkX graph
        :return: None
        """
        g = nx.Graph()
        # For each vertex and edge in molecule graph add node and edge in NetworkX graph
        for n in self.structure.vertices():
            g.add_node(self.structure.position_of_vertex(n), element=n.element)
        for e in self.structure.edges():
            if isinstance(e, m.Bond):
                if e.single:
                    g.add_edge(self.structure.endpoints_position(e)[0], self.structure.endpoints_position(e)[1], type='single')
                elif e.double:
                    g.add_edge(self.structure.endpoints_position(e)[0], self.structure.endpoints_position(e)[1], type='double')
                elif e.triple:
                    g.add_edge(self.structure.endpoints_position(e)[0], self.structure.endpoints_position(e)[1], type='triple')
                elif e.quadruple:
                    g.add_edge(self.structure.endpoints_position(e)[0], self.structure.endpoints_position(e)[1], type='quadruple')
                elif e.aromatic:
                    g.add_edge(self.structure.endpoints_position(e)[0], self.structure.endpoints_position(e)[1], type='aromatic')
        self.nx_graph = g

if __name__ == '__main__':
    struct = PathStructure(smiles_parser.Parser().parse_smiles('C1CCC1'))
    struct.create_nx_graph()