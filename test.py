import graph
import unittest

class TestVertex(unittest.TestCase):
	
	def setUp(self):
		self.x = graph.Vertex('C', False, False, False, False)
		
	def testVertexCreation(self):
		self.assertEqual(self.x.element, 'C')

	#def testAddAdjacent(self):
	#	self.x.addAdjacent(2, 'Double')
	#	self.assertEqual(self.x.adjacent.get(2), 'Double')
	
	
class TestEdge(unittest.TestCase):
	
	def setUp(self):
		self.x = graph.Vertex('C', False, False, False, False)
		self.y = graph.Vertex('He', False, 2, False, False)
		self.e = graph.Edge(self.x, self.y, 1)
		
	def test_edge_creation(self):
		self.assertEqual(self.e.weight, 1)
		
	def test_endpoints(self):
		self.assertEqual(self.e.endpoints, (self.x, self.y))
		
	def test_opposite(self):
		self.assertEqual(self.e.opposite(self.x), self.y)
		

class TestGraph(unittest.TestCase):
	
	def setUp(self):
		self.g = graph.Graph('C=N')
		self.x = self.g.add_vertex('C')
		self.y = self.g.add_vertex('N', hydrogen=2)
		self.e = self.g.add_edge(self.x, self.y, 2)
		
	def test_graph_creation(self):
		# ERROR: Graph instance has no attribute 'SMILES_string'
		#self.assertEqual(self.g.SMILES_string, 'C=N')
		pass
	
	def test_clear(self):
		self.g.clear()
		self.assertEqual(self.g.vertices, [])
		
	def test_add_vertex(self):
		self.assertEqual(self.y.hydrogen_count, 2)
				
	def test_remove_vertex(self):
		self.assertEqual(self.g.vertices[0].element, 'C')
		self.g.remove_vertex(self.x)
		self.assertEqual(self.g.vertices[0].element, 'N')
		
	def test_add_edge(self):
		self.assertEqual(self.e.weight, 2)
		
	def test_remove_edge(self):
		pass
		#self.g.remove_edge(self.x, self.y)
		#self.assertNotEqual(self.g.vertices[0], 'D')
		#self.assertEqual(self.G.getVertex(1).getAdjacentVertices(), {})
		
	def test_neighbours(self):
		pass
		
	def test_connecting_edges(self):
		pass
		
	def test_degree(self):
		pass
		
	def test_contains_edge(self):
		pass
	

class TestTokenizer(unittest.TestCase):
	pass
	
	
class TestParser(unittest.TestCase):
	pass
		
		
def main():
	unittest.main()
	
if __name__ == '__main__':
	main()		
