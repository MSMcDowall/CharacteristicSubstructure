import molecule
import igraph

def draw_molecule(mole):
    edge_list = []
    for e in mole.edges:
        if e is molecule.SingleBond:
            edge_list.append(e.endpoints_position())
        elif e is molecule.DoubleBond:
            edge_list.append(e.endpoints_position() * 2)
        elif e is molecule.TripleBond:
            edge_list.append(e.endpoints_position() * 3)
        elif e is molecule.QuadrupleBond:
            edge_list.append(e.endpoints_position() * 4)
    g = igraph.Graph()
    g.add_vertices(mole.size)
    g.add_edges(edge_list)
    vertex_labels = []
    for v in mole.vertices:
        vertex_labels.append(str(v.element))
    g.vs["label"] = vertex_labels
    return g



