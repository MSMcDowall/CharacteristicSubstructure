class Vertex:
	def __init__(self, id, type):
		self.id = id
		self.type = type
		self.adjacent = {}

	def getId(self):
		return self.id

	def getType(self):
		return self.type

	def getAdjacentVertices(self):
		return self.adjacent.keys()

	def addAdjacent(self, otherId, weight):
		self.adjacent[otherId] = weight
		
		
