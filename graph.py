# A single Vertex of Graph which includes its element
class Vertex(object):
    def __init__(self, element):
        self._element = element

    # Return the element of the vertex
    @property
    def element(self):
        return self._element

    # Create a hash of the vertex object so it can be used as a key
    def __hash__(self):
        return hash(id(self))


# A single Edge of Graph which includes the weight of the chemical element
class Edge(object):
    def __init__(self, origin, destination, element=None):
        self._origin = origin
        self._destination = destination
        self._element = element

    # Return the endpoints of the edge as a tuple
    @property
    def endpoints(self):
        return (self._origin, self._destination)

    # Return the element of the edge
    @property
    def element(self):
        return self._element

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
class Graph(object):
    def __init__(self):
        # Dictionary of vertices that have been added to the Graph which maps them to their adjacent vertices
        self._vertices = {}
        # The number of vertices within the graph
        self._vertex_count = 0

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

    # Clears the graph
    def clear(self):
        self._vertices = {}
        self._vertex_count = 0

    # Creates a new vertex object and assigns it a dictionary which will contain all adjacent vertices and edges
    # Method is split in two to aid in the inheritance of this method
    def add_vertex(self, element):
        new_vertex = Vertex(element)
        self.vertex_to_graph(new_vertex)
        return new_vertex

    def vertex_to_graph(self, vertex):
        self._vertices[vertex] = {}
        self._vertex_count = + 1

    # Deletes the vertex object
    def remove_vertex(self, vertex):
        if vertex in self._vertices:
            # Delete the vertex object from the dictionaries of vertices which are adjacent to it
            # del(self._vertices[adjacent].get(vertex) for adjacent in self._vertices[vertex])
            del self._vertices[vertex]

    # Add a weighted edge between the vertices at the two given positions
    # Method is split in two to aid in the inheritance of this method
    def add_edge(self, first_vertex, second_vertex, element=None):
        new_edge = Edge(first_vertex, second_vertex, element)
        self.edge_to_graph(first_vertex, second_vertex, new_edge)
        return new_edge

    def edge_to_graph(self, first_vertex, second_vertex, edge):
        self._vertices[first_vertex][second_vertex] = edge
        self._vertices[second_vertex][first_vertex] = edge

    # Remove the edge between the V=vertices at the two given positions
    def remove_edge(self, first_vertex, second_vertex):
        del self._vertices[first_vertex][second_vertex]
        del self._vertices[second_vertex][first_vertex]


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
