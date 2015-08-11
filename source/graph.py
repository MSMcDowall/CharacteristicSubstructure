import copy
#from draw_molecule import draw_molecule as draw

# A single Vertex of Graph which includes its element
class Vertex(object):
    def __init__(self, element):
        self.element = element
        self.position = 0       # The position of the vertex in the graph
        self.visited = False    # An indicator that is used in DFS

    # Create a hash of the vertex object so it can be used as a key
    def __hash__(self):
        return hash(id(self))


# A single Edge of Graph which includes the weight of the chemical element
class Edge(object):
    def __init__(self, origin, destination, element=None):
        self._origin = origin
        self._destination = destination
        self.element = element

    # Return the endpoints of the edge as a tuple
    @property
    def endpoints(self):
        return (self._origin, self._destination)

    # Return the position of the two endpoints
    def endpoints_position(self):
        return (self._origin.position, self._destination.position)

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

    def __str__(self):
        return 'join %s and %s' % (self._origin, self._destination)


# A Graph object which maintains a dictionary of Vertex objects		
class Graph(object):
    def __init__(self):
        self.adjacency_dictionary = {}         # Dictionary of vertices and maps them to their adjacent vertices
        self.size = 0               # The number of vertices within the graph
        self.paths = []             # A collection of all the path strings as well as the vertices that are involved

    # Returns all the vertices in the graph
    def vertices(self):
        return self.adjacency_dictionary.keys()

    # Returns all the edges of the graph
    def edges(self):
        edges = set()
        for adjacentVertex in self.adjacency_dictionary.values():
            edges.update(adjacentVertex.values())
        return edges

    # Clears the graph
    def clear(self):
        self.adjacency_dictionary = {}
        self.size = 0

    # Creates a new vertex object and assigns it a dictionary which will contain all adjacent vertices and edges
    # Method is split in two to aid in the inheritance of this method
    def add_vertex(self, element):
        new_vertex = Vertex(element)
        self.vertex_to_graph(new_vertex)
        return new_vertex

    # Adds a vertex that has already been created to the graph 
    def vertex_to_graph(self, vertex):
        self.adjacency_dictionary[vertex] = {}
        vertex.position = self.size
        self.size += 1

    # Deletes the vertex object
    def remove_vertex(self, vertex):
        # Delete the vertex object from the dictionaries of vertices which are adjacent to it
        for neighbour in self.neighbours(vertex):
            del self.adjacency_dictionary[neighbour][vertex]
        del self.adjacency_dictionary[vertex]

    # Swap out a vertex and replace it with the vertex object provided
    def swap_vertex(self, old_vertex, new_vertex):
        self.adjacency_dictionary[new_vertex] = copy.copy(self.adjacency_dictionary[old_vertex])
        print 'in swapping'
        print self.adjacency_dictionary
        # draw(self)
        print 'neighbours of next'
        print repr(old_vertex)
        print self.neighbours(old_vertex)
        for neighbour in self.neighbours(old_vertex):
            print 'neighbour'
            print repr(neighbour)
            print self.neighbours(neighbour)
            self.adjacency_dictionary[neighbour][new_vertex] = copy.copy(self.adjacency_dictionary[neighbour][old_vertex])

        # Change any instances where the old vertex appears in the list of paths
        if self.paths:
            path_copy = copy.copy(self.paths)
            for path_tuple in self.paths:
                if old_vertex in path_tuple[1]:
                    tuple_index = self.paths.index(path_tuple)
                    vertex_index = path_tuple[1].index(old_vertex)
                    vertices_list = copy.copy(path_tuple[1])
                    vertices_list[vertex_index] = new_vertex
                    path_copy[tuple_index] = (path_tuple[0], vertices_list)
            self.paths = copy.copy(path_copy)
        self.remove_vertex(old_vertex)

    # Add a weighted edge between the vertices at the two given positions
    # Method is split in two to aid in the inheritance of this method
    def add_edge(self, first_vertex, second_vertex, element=None):
        new_edge = Edge(first_vertex, second_vertex, element)
        self.edge_to_graph(first_vertex, second_vertex, new_edge)
        return new_edge

    # Adds a bond that has already been created to the graph
    def edge_to_graph(self, first_vertex, second_vertex, edge):
        self.adjacency_dictionary[first_vertex][second_vertex] = edge
        self.adjacency_dictionary[second_vertex][first_vertex] = edge

    # Remove the edge between the V=vertices at the two given positions
    def remove_edge(self, first_vertex, second_vertex):
        del self.adjacency_dictionary[first_vertex][second_vertex]
        del self.adjacency_dictionary[second_vertex][first_vertex]

    # # Create a printable string version of each vertex dictionary
    # def dictionary_string(self, vertex):
    #     dictionary_string = {}
    #     for key, value in self.adjacency_dictionary[vertex].iteritems():
    #         dictionary_string[str(key)] = str(value)
    #     return dictionary_string

    # Returns all the vertices which are adjacent to the vertex
    def neighbours(self, vertex):
        return self.adjacency_dictionary[vertex].keys()

    # Returns all the edges which connect to the vertex
    def connecting_edges(self, vertex):
        return self.adjacency_dictionary[vertex].values()

    # Returns the number of adjacent vertices to the given vertex
    def degree(self, vertex):
        return len(self.adjacency_dictionary[vertex])

    # Tests if there is an edge between the two vertices
    def contains_edge(self, first_vertex, second_vertex):
        if second_vertex in self.adjacency_dictionary[first_vertex]:
            return self.adjacency_dictionary[first_vertex][second_vertex]
        else:
            return False

    # A depth first search which visits each node in turn
    # A set is then compiled of all the paths within the graph
    # Algorithm structure from Handbook of Graph Theory, Gross & Yellen
    def find_all_paths(self):
        completed = []      # The nodes which have acted as a root for the search
        all_paths = {}      # The dictionary of all the paths and their lengths (dict removes duplicate paths)
        for v in self.vertices():
            if v not in completed:
                for w in self.vertices():
                    w.visited = False       # Start search anew for each root
                path_stack = []             # Used to create the string of the path
                position_stack = []         # Used to create the string of positions
                self.find(v, path_stack, position_stack, all_paths)
                completed.append(v)
        return all_paths

    # The recursive element of the depth first search
    def find(self, v, path_stack, position_stack, all_paths):
        v.visited = True
        path_stack.append(v.element + '-')
        position_stack.append(v)
        for w in self.adjacency_dictionary[v].keys():
            if not w.visited:
                self.find(w, path_stack, position_stack, all_paths)
        letters = ''.join(path_stack)
        path = letters[:-1]    # Remove final dash
        positions = list(position_stack)
        all_paths[path] = len(letters)/2
        self.paths.append((path, positions))
        path_stack.pop()
        position_stack.pop()

if __name__ == '__main__':
    G = Graph()
    a = G.add_vertex('C')
    b = G.add_vertex('N')
    c = G.add_vertex('C')
    G.add_edge(a, b)
    G.add_edge(b, c)
    print G.adjacency_dictionary
    new = Vertex('new')
    print new
    G.swap_vertex(b, new)
    print 'vertices'
    for v in G.adjacency_dictionary:
        print v
    print G.neighbours(new)
    print G.neighbours(a)