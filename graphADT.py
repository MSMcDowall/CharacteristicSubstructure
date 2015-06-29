# A single Vertex of Graph which includes its chemical element and a dictionary of adjacent vertices
class Vertex:
    def __init__(self, element, isotope, hydrogen, charge, aromatic):
        self._element = element
        self._isotope = isotope
        self._hydrogen = hydrogen
        self._charge = charge
        self._aromatic = aromatic

    # Return the element coloring of the vertex
    @property
    def element(self):
        return self._element

    # Return the isotope number of the vertex
    @property
    def isotope(self):
        return self._isotope

    # Return the number of hydrogens attached to the vertex
    @property
    def hydrogen_count(self):
        return self._hydrogen

    # Return the charge of the element
    @property
    def charge(self):
        return self._charge

    # Return if the element is aromatic
    @property
    def artomatic(self):
        return self._aromatic

    # Create a hash of the vertex object so it can be used as a key
    def __hash__(self):
        return hash(id(self))


# A single Edge of Graph which includes the weight of the chemical element
class Edge:
    def __init__(self, origin, destination, weight):
        self._origin = origin
        self._destination = destination
        self._weight = weight

    # Return the weight of the an edge
    @property
    def weight(self):
        return self._weight

    # Return the endpoints of the edge as a tuple
    @property
    def endpoints(self):
        return (self._origin, self._destination)

    # Return the opposite endpoint of the edge to the one given as a parameter
    def opposite(self, vertex):
        if vertex == self._origin:
            return self._destination
        if vertex == self._destination:
            return self._origin
        else:
            return None

    # Create a hash of the edge object
    def __hash__(self):
        return hash((self._origin, self._destination))


# A Graph object which maintains a dictionary of Vertex objects		
class Graph:
    def __init__(self, string):
        # Dictionary of vertices that have been added to the Graph which maps them to their adjacent vertices
        self._vertices = {}
        # Stores the elements whose bonds have been broken within a ring
        # The break signifying number is used as a key
        self._breakpoints = {}
        # The number of vertices within the graph
        self._vertex_count = 0
        # The original SMILES string
        self._SMILES_string = string

    # Returns all the vertices in the graph
    @property
    def vertices(self):
        return self._vertices.keys()

    # Returns all the edges of the graph
    @property
    def edges(self):
        edges = set()
        for adjacentVertex in self._vertices.values():
            edges.update(adjacentVertex.values())
        return edges

    # Returns the number of vertices in the graph
    @property
    def size(self):
        return len(self._vertices)
        
	# Return the original SMILES string
	@property
	def SMILES_string(self):
		return self._SMILES_string

    # Clears the graph
    def clear(self):
        self._vertices = {}
        self._breakpoints = {}
        self._vertex_count = 0

    # Creates a new vertex object and assigns it a dictionary which will contain all adjacent vertices and edges
    def add_vertex(self, element, isotope=False, hydrogen=False, charge=False, aromatic=False):
        new_vertex = Vertex(element, isotope, hydrogen, charge, aromatic)
        self._vertices[new_vertex] = {}
        self._vertex_count = + 1
        return new_vertex

    # Deletes the vertex object
    def remove_vertex(self, vertex):
        if vertex in self._vertices:
            # Delete the vertex object from the dictionaries of vertices which are adjacent to it
            # del(self._vertices[adjacent].get(vertex) for adjacent in self._vertices[vertex])
            del self._vertices[vertex]

    # Add a weighted edge between the vertices at the two given positions
    def add_edge(self, first_vertex, second_vertex, weight=1):
        if first_vertex in self._vertices and second_vertex in self._vertices:
            newEdge = Edge(first_vertex, second_vertex, weight)
            self._vertices[first_vertex][second_vertex] = newEdge
            self._vertices[second_vertex][first_vertex] = newEdge
            return newEdge
        else:
            print 'Vertex missing'

    # Remove the edge between the V=vertices at the two given positions
    def remove_edge(self, first_vertex, second_vertex):
        if first_vertex in self._vertices and second_vertex in self._vertices:
            del self._vertices[first_vertex][second_vertex]
            del self._vertices[second_vertex][first_vertex]
        else:
            print 'Vertex missing'

    # Returns all the vertices which are adjacent to the vertex
    def neighbours(self, vertex):
        return self._vertices[vertex].keys()

    # Returns all the edges which connect to the vertex
    def connecting_edges(self, vertex):
        return self._vertices[vertex].values()

    # Returns the number of adjacent vertices to the given vertex
    def degree(self, vertex):
        return len(self._vertices[vertex])

    # Tests if there is an edge between the two vertices
    def contains_edge(self, first_vertex, second_vertex):
        if second_vertex in self._vertices[first_vertex]:
            return True
        else:
            return False
