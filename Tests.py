import GraphADT
import unittest

class TestVertex(unittest.TestCase):
	
	def setUp(self):
		self.x = GraphADT.Vertex(1, 'C')
		
	def testVertexCreation(self):
		self.assertEqual(self.x.getPosition(), 1)

	def testAddAdjacent(self):
		self.x.addAdjacent(2, 'Double')
		self.assertEqual(self.x.adjacent.get(2), 'Double')
	
		
class TestGraph(unittest.TestCase):
	
	def setUp(self):
		self.G = GraphADT.Graph()
		self.G.addVertex('Carbon', 0)
		self.G.addVertex('Nitrogen', 1)
		self.G.addEdge(0, 1, 'D')
		
	def testAddVertex(self):
		self.assertEqual(self.G.getVertex(0).getElement(), 'Carbon')
		self.assertEqual(self.G.getVertex(1).getElement(), 'Nitrogen')
	
	def testRemoveVertex(self):
		self.G.removeVertex(0)
		self.assertIsNone(self.G.getVertex(0))
		
	def testAddEdge(self):
		
		self.assertEqual(self.G.getVertex(0).getAdjacentVertices(), {1: 'D'})
		self.assertNotEqual(self.G.getVertex(1).getAdjacentVertices(), {})
		
	def testRemoveEdge(self):
		self.G.removeEdge(1, 0)
		self.assertNotEqual(self.G.getVertex(0).getAdjacentVertices(), {1: 'D'})
		self.assertEqual(self.G.getVertex(1).getAdjacentVertices(), {})
				
	
		
def main():
	unittest.main()
	
if __name__ == '__main__':
	main()		
