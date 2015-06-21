class Vertex:
	def __init__(self, position, element):
		self.position = position
		self.element = element
		self.adjacent = {}

	def getPosition(self):
		return self.position

	def getElement(self):
		return self.element
		
	def addAdjacent(self, otherPosition, weight):
		self.adjacent[otherPosition] = weight

	def getAdjacentVertices(self):
		return self.adjacent.keys()

		
class Graph:
	def __init__(self, name):
		self.name = name
		self.vertices = []
		self.breakpoints = {}
		self.numVertices = 0
		
	def addVertex(self, element):
		self.vertices.append(Vertex(self.numVertices, element))
		self.numVertices =+ 1
	
	def getVertex(self, index):
		return self.vertices[index]
		
				
