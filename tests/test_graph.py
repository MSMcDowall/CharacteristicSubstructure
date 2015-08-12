import source.graph as graph
import unittest


class VertexTestCase(unittest.TestCase):
    def setUp(self):
        self.x = graph.Vertex('C')

    def test_vertex_creation(self):
        self.assertEqual(self.x.element, 'C')


class EdgeTestCase(unittest.TestCase):
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


class GraphTestCase(unittest.TestCase):
    def setUp(self):
        self.g = graph.Graph()
        self.x = self.g.add_vertex('C')
        self.y = self.g.add_vertex('N')
        self.e = self.g.add_edge(self.x, self.y, 2)

    def test_graph_creation(self):
        pass

    def test_clear(self):
        self.g.clear()
        self.assertEqual(self.g.vertices(), [])

    def test_add_vertex(self):
        self.assertEqual(self.y.element, 'N')

    def test_remove_vertex(self):
        self.assertEqual(self.g.vertices()[0].element, 'C')
        self.g.remove_vertex(self.x)
        self.assertEqual(self.g.vertices()[0].element, 'N')

    def test_add_edge(self):
        self.assertEqual(self.e.element, 2)

    def test_remove_edge(self):
        self.assertEqual(self.g.adjacency_dictionary[self.x][self.y].element, 2)
        self.g.remove_edge(self.x, self.y)
        self.assertEqual(self.g.adjacency_dictionary[self.x], {})

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

    def test_swap_vertex(self):
        self.a = graph.Vertex('A')
        self.b = self.g.add_vertex('B')
        self.f = self.g.add_edge(self.x, self.b, 'f')
        self.assertIn(self.b, self.g.adjacency_dictionary[self.x])
        self.g.swap_vertex(self.b, self.a)
        self.assertIn(self.a, self.g.adjacency_dictionary[self.x])
        self.assertIn(self.x, self.g.adjacency_dictionary[self.a])
        self.assertEqual(self.g.adjacency_dictionary[self.a][self.x].element, 'f')
        self.assertNotIn(self.b, self.g.adjacency_dictionary)

    def test_swap_path_vertices(self):
        self.a = graph.Vertex('C')
        self.g.find_all_paths()
        self.assertIn(('C', [self.x]), self.g.paths)
        self.assertIn(('C-N', [self.x, self.y]), self.g.paths)
        self.g.swap_vertex(self.x, self.a)
        self.assertIn(('C', [self.a]), self.g.paths)
        self.assertIn(('C-N', [self.a, self.y]), self.g.paths)
        self.assertNotIn(('C', [self.x]), self.g.paths)
        self.assertNotIn(('C-N', [self.x, self.y]), self.g.paths)

    def find_all_paths(self):
        self.g.find_all_paths()
        self.assertIn(('C', [self.x]), self.g.paths)
        self.assertIn(('N', [self.y]), self.g.paths)
        self.assertIn(('C-N', [self.x, self.y]), self.g.paths)

if __name__ == '__main__':
    unittest.main()
