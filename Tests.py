import GraphADT
import unittest

class TestVertex(unittest.TestCase):
	
	def setUp(self):
		self.x = GraphADT.Vertex(1, 'C')
		
	def testVertexCreation(self):
		self.assertEqual(self.x.getPosition(), 1)

	def testAddAdjacent(self):
		self.x.addAdjacent(2, 'Double')
		self.assertEqual(str(self.x.adjacent.get(2)), 'Double')
	
		
class TestGraph(unittest.TestCase):
	
	def setUp(self):
		self.G = GraphADT.Graph('Bob')
		
	def testAddVertex(self):
		self.G.addVertex('Carbon')
		self.assertEqual(self.G.getVertex(0).getElement(), 'Carbon')
		self.G.addVertex('Nitrogen')
		self.assertEqual(self.G.getVertex(1).getElement(), 'Nitrogen')
		
def main():
	unittest.main()
	
if __name__ == '__main__':
	main()		
