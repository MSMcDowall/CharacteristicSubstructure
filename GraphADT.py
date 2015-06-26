# A single Vertex of Graph which includes its chemical element and a dictionary of adjacent vertices
class Vertex:
	def __init__(self, position, element):
		self.position = position
		self.element = element
		self.adjacent = {}

	# The postition of the Vertex in its Graph is used as an identifier
	def getPosition(self):
		return self.position

	# Return the element coloring of the vertex
	def getElement(self):
		return self.element

	# Return the positions of each of the adjacent vertices in the Graph
	def getAdjacentVertices(self):
		return self.adjacent
		
# A single Edge of Graph which includes the weight of the chemical element
class Edge:
	def __init__(self, origin, destination, weight):
		self.origin = origin
		self.destination = destination
		self.weight = weight
		
	

# A Graph object which maintains a dictionary of Vertex objects		
class Graph:
	def __init__(self):
		# Dictionary of vertices that have been added to the Graph
		self.vertices = {}
		# Stores the elements whose bonds have been broken within a ring
		# The break signifying number is used as a key
		self.breakpoints = {}
		# The number of vertices within the graph
		self.numVertices = 0
	
	# Creates a new Vertex object and adds it to the list at the given index postion	
	def addVertex(self, element, position):
		newVertex = Vertex(position, element)
		self.vertices[position] = newVertex
		self.numVertices =+ 1
	
	# Returns the Vertex object which is at the given position if it exists
	def getVertex(self, position):
		if position in self.vertices:
			return self.vertices[position]
		else:
			return None
	
	# Deletes the Vertex object at the given position
	def removeVertex(self, position):
		if position in self.vertices:
			del self.vertices[position]
			
	# Add a weighted Edge between the Vertices at the two given positions
	def addEdge(self, firstPosition, secondPosition, weight='S'):
		if firstPosition in self.vertices and secondPosition in self.vertices:
			self.vertices[firstPosition].addAdjacent(secondPosition, weight)
			self.vertices[secondPosition].addAdjacent(firstPosition, weight)
		else: 
			print 'Vertex missing'
		
	# Remove the edge between the Vertices at the two given positions
	def removeEdge(self, firstPosition, secondPosition):
		if firstPosition in self.vertices and secondPosition in self.vertices:
			del self.vertices[firstPosition].adjacent[secondPosition]
			del self.vertices[secondPosition].adjacent[firstPosition]
		else:
			print 'Vertex missing'
