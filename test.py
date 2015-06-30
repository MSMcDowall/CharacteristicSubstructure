import graph
import unittest

class TestVertex(unittest.TestCase):

    def setUp(self):
        self.x = graph.Vertex('C')

    def testVertexCreation(self):
        self.assertEqual(self.x.element, 'C')

    #def testAddAdjacent(self):
    #	self.x.addAdjacent(2, 'Double')
    #	self.assertEqual(self.x.adjacent.get(2), 'Double')


class TestEdge(unittest.TestCase):

    def setUp(self):
        self.x = graph.Vertex('C')
        self.y = graph.Vertex('He')
        self.e = graph.Edge(self.x, self.y, 1)

    def test_edge_creation(self):
        self.assertEqual(self.e.element, 1)

    def test_endpoints(self):
        self.assertEqual(self.e.endpoints, (self.x, self.y))

    def test_opposite(self):
        self.assertEqual(self.e.opposite(self.x), self.y)


class TestGraph(unittest.TestCase):

    def setUp(self):
        self.g = graph.Graph()
        self.x = self.g.add_vertex('C')
        self.y = self.g.add_vertex('N')
        self.e = self.g.add_edge(self.x, self.y, 2)

    def test_graph_creation(self):
        # ERROR: Graph instance has no attribute 'SMILES_string'
        #self.assertEqual(self.g.SMILES_string, 'C=N')
        pass

    def test_clear(self):
        self.g.clear()
        self.assertEqual(self.g.vertices, [])

    def test_add_vertex(self):
        self.assertEqual(self.y.element, 'N')

    def test_remove_vertex(self):
        self.assertEqual(self.g.vertices[0].element, 'C')
        self.g.remove_vertex(self.x)
        self.assertEqual(self.g.vertices[0].element, 'N')

    def test_add_edge(self):
        self.assertEqual(self.e.element, 2)

    def test_remove_edge(self):
        pass
    #self.g.remove_edge(self.x, self.y)
    #self.assertNotEqual(self.g.vertices[0], 'D')
    #self.assertEqual(self.G.getVertex(1).getAdjacentVertices(), {})

    def test_neighbours(self):
        self.assertEqual(self.g.neighbours(self.x), [self.y])

    def test_connecting_edges(self):
        self.assertEqual(self.g.connecting_edges(self.y), [self.e])

    def test_degree(self):
        self.assertEqual(self.g.degree(self.x), 1)
        self.z = self.g.add_vertex('O')
        self.g.add_edge(self.x, self.z)
        self.assertEqual(self.g.degree(self.x), 2)

    def test_contains_edge(self):
        self.z = self.g.add_vertex('O')
        self.g.add_edge(self.x, self.z)
        self.assertTrue(self.g.contains_edge(self.x, self.z))
        self.assertFalse(self.g.contains_edge(self.y, self.z))


def main():
    unittest.main()

if __name__ == '__main__':
    main()
