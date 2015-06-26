import graphADT
import unittest

class TestVertex(unittest.TestCase):
	
	def setUp(self):
		self.x = graphADT.Vertex(1, 'C')
		
	def testVertexCreation(self):
		self.assertEqual(self.x.getPosition(), 1)

	def testAddAdjacent(self):
		self.x.addAdjacent(2, 'Double')
		self.assertEqual(self.x.adjacent.get(2), 'Double')
	
		
class TestGraph(unittest.TestCase):
	
	def setUp(self):
		self.G = graphADT.Graph()
		self.G.add_vertex('Carbon', 0)
		self.G.add_vertex('Nitrogen', 1)
		self.G.add_edge(0, 1, 'D')
		
	def testAddVertex(self):
		self.assertEqual(self.G.getVertex(0).element(), 'Carbon')
		self.assertEqual(self.G.getVertex(1).element(), 'Nitrogen')
	
	def testRemoveVertex(self):
		self.G.remove_vertex(0)
		self.assertIsNone(self.G.getVertex(0))
		
	def testAddEdge(self):
		
		self.assertEqual(self.G.getVertex(0).getAdjacentVertices(), {1: 'D'})
		self.assertNotEqual(self.G.getVertex(1).getAdjacentVertices(), {})
		
	def testRemoveEdge(self):
		self.G.remove_edge(1, 0)
		self.assertNotEqual(self.G.getVertex(0).getAdjacentVertices(), {1: 'D'})
		self.assertEqual(self.G.getVertex(1).getAdjacentVertices(), {})
		
class TestEdge(unittest.TestCase):
	
	pass
				
	
		
def main():
	unittest.main()
	
if __name__ == '__main__':
	main()		
