# A single Vertex of Graph which includes its chemical element and a dictionary of adjacent vertices
class Vertex:
	def __init__(self, element, isotope, hydrogen, charge):
		self._element = element
		self._isotope = isotope
		self._hydrogen = hydrogen
		self._charge = charge

	# Return the element coloring of the vertex
	def getElement(self):
		return self._element

	# Return the isotope number of the vertex
	def getIsotope(self):
		return self._isotope
		
	# Return the number of hydrogens attached to the vertex
	def getHydrogenCount(self):
		return self._hydrogen
		
	# Return the charge of the element
	def getCharge(self):
		return self._charge
		
	# Create a hash of the vertex object so it can be used as a key
	def __hash__(self):
		return hash(id(self))
		
# A single Edge of Graph which includes the weight of the chemical element
class Edge:
	def __init__(self, origin, destination, weight):
		self._origin = origin
		self._destination = destination
		self._weight = weight
	
	# Return the weight of the 
	def getWeight(self):
		return self._weight
	
	# Return the endpoints of the edge as a tuple	
	def endpoints(self):
		return (self._origin, self._destination)
		
	# Return the opposite endpoint of the edge to the one given as a parameter
	def opposite(self, vertex):
		if vertex == self._origin:
			return self._destination
		if vertex == self._destination;
			return self._origin
		else:
			return null
			
	# Create a hash of the edge object
	def __hash__(self):
		return hash((self._origin, self._destination))
		

# A Graph object which maintains a dictionary of Vertex objects		
class Graph:
	def __init__(self):
		# Dictionary of vertices that have been added to the Graph which maps them to their adjacent vertices
		self._vertices = {}
		# Stores the elements whose bonds have been broken within a ring
		# The break signifying number is used as a key
		self._breakpoints = {}
		# The number of vertices within the graph
		self._numVertices = 0
	
	# Returns all the vertices in the graph
	def getVertices(self):
		return self._vertices.keys()
	
	# Returns all the edges of the graph	
	def getEdges(self):
		edges = set()
		for adjacentVertex in self._vertices.values():
			edges.update(adjacentVertex.values())
		return edges
		
		
	# Returns the number of vertices in the graph
	def size(self):
		return len(self._vertices)
		
	# Clears the graph
	def clear(self):
		self = Graph()
	
	# Creates a new Vertex object and adds it to the list at the given index postion	
	def addVertex(self, element, isotope='natural', hydrogen='lowestValence', charge='neutral'):
		newVertex = Vertex(element, isotope, hydrogen, charge)
		self.vertices[newVertex] = {}
		self.numVertices =+ 1
		return newVertex
	
	# Deletes the Vertex object
	def removeVertex(self, vertex):
		if vertex in self._vertices:
			# Delete the vertex object from the dictionaries of vertices which are adjacent to it
			del(self._vertices[adjacent].get(vertex) for adjacent in self._vertices[vertex])
			del self._vertices[vertex]
			
			
	# Add a weighted Edge between the Vertices at the two given positions
	def addEdge(self, firstPosition, secondPosition, weight='S'):
		if firstPosition in self._vertices and secondPosition in self._vertices:
			newEdge = self.Edge(firstPosition, secondPosition, weight)
			self._vertices[firstPosition][secondPostion] = newEdge
			self._vertices[secondPosition][firstPosition] = newEdge
			return newEdge
		else: 
			print 'Vertex missing'
		
	# Remove the edge between the Vertices at the two given positions
	def removeEdge(self, firstPosition, secondPosition):
		if firstPosition in self.vertices and secondPosition in self.vertices:
			del self._vertices[firstPosition][secondPosition]
			del self._vertices[secondPosition][firstPosition]
		else:
			print 'Vertex missing'
				
	#neighbours(vertex)
	def neighbours(self, vertex):
		return self._vertices[vertex].keys()
		
	#connectingEdges(vertex)
	def connectingEdges(self, vertex):
		return self._vertices[vertex].values()
		
	#degree(vertex)
	def degree(self, vertex):
		return len(self._vertices[vertex])
	
	#containsEdge(vertex1, vertex2)
	def contains_edge(self, first_vertex, second_vertex):
		if second_vertex in self._vertices[first_vertex]:
			return True
		else:
			return False
