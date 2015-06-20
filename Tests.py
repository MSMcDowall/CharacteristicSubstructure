from GraphADT import Vertex
import unittest

class TestGraphADt(unittest.TestCase):
	
	def setUp(self):
		self.x = Vertex(1, 'C')
		self.x.addAdjacent(2, 'Double')
		
	def testVertexCreation(self):
		self.assertEqual(self.x.getId(), 1)

	def testAddAdjacent(self):
		self.assertEqual(str(self.x.adjacent.get(2)), 'Double')
		
def main():
	unittest.main()
	
if __name__ == '__main__':
	main()		
