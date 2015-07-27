import molecule
import parser
import networkx as nx
import matplotlib.pyplot as plt


def draw_molecule(mole):
    # Create a new NetworkX graph
    g = nx.Graph()

    # For each vertex and edge in molecule graph add node and edge in NetworkX graph
    for n in mole.vertices:
        g.add_node(n.position, element=n.element)
    for e in mole.edges:
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

    # Set the layout
    pos = nx.spring_layout(g, iterations=20)
    # Display the element type and edge type as labels
    labels = dict((n,d['element']) for n,d in g.nodes(data=True))
    edge_labels = dict(((u,v),d['type']) for u,v,d in g.edges(data=True))
    # Add the labels to the graph
    nx.draw(g, pos=pos, node_color='w')
    nx.draw_networkx_labels(g, pos=pos, labels=labels)
    nx.draw_networkx_edge_labels(g, pos=pos, edge_labels=edge_labels)
    # Display the completed graph
    plt.show()
    return g

if __name__ == '__main__':
    mole = parser.Parser().parse_smiles('ClBrCB')
    print mole.vertices
    graph = draw_molecule(mole)


