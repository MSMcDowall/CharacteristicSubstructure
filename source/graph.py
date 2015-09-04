from copy import copy
#from draw_molecule import draw_molecule as draw

"""
This module contains a graph implemented with an adjacency dictionary
and the edge and vertex objects which can be added to the graph.
Based on the Graph ADT in 'Algorithms and Data Structures in Python'
by  Michael T. Goodrich, Roberto Tamassia, Michael H. Goldwasser
"""


class Vertex(object):
    """
    A single vertex of a graph which includes its element.

    The vertex does not hold information about its connections; this information is stored in the graph
    :param element: a string which is the element of the vertex
    position: an integer which is used to signify the order in which the vertex was added to the graph
            it can be used as a simple unique identifier for the vertex within the graph
    visited: a boolean flag which is used in the find_all_paths method of a graph
    """
    def __init__(self, element):
        self.element = element
        self.visited = False

    def __hash__(self):
        return hash(id(self))


class Edge(object):
    """
    A single Edge of a Graph which includes the two vertices that it connects and its element.

    :param origin: one of the vertex objects which the edge connects
    :param destination: the other vertex which the edge connects
    :param element: a string which is the element of the edge
    """
    def __init__(self, origin, destination, element=None):
        self._origin = origin
        self._destination = destination
        self.element = element

    @property
    def endpoints(self):
        """
        Gives the endpoints of the edge
        :return: a tuple of the two endpoints
        """
        return self._origin, self._destination

    def opposite(self, vertex):
        """
        Gives the vertex object which lies opposite to the given vertex on the edge
        :param vertex: the vertex object whose opposite is to be found
        :return: the vertex object that is the other endpoint in the edge
        """
        if vertex == self._origin:
            return self._destination
        if vertex == self._destination:
            return self._origin
        else:
            return None

    def __hash__(self):
        return hash((self._origin, self._destination))

    def __str__(self):
        return 'join %s and %s' % (self._origin, self._destination)


class Graph(object):
    """
    A Graph which contains an adjacency dictionary that contains the information about its vertices and edges.

    adjacency_dictionary: a dictionary that takes the form of {vertex: {adjacent vertex: connecting edge}}
    size: an integer gives the number of individual vertices in the graph
    paths: a list of all of the paths which are contained in the graph in tuples with vertices that are in the path
    """
    def __init__(self):
        self.adjacency_dictionary = {}
        self.size = 0
        self.paths = []
        # list of vertices, index in list acts as position indicator for isomorphism checking
        self.positions = []

    def vertices(self):
        """
        Returns all the vertices of the graph

        Performance: O(n)
        :return: a list of all the vertex objects present in the graph
        """
        return self.adjacency_dictionary.keys()

    def edges(self):
        """
        Returns all the edges of the graph in a set to remove duplicates

        Performance: O(m)
        :return: a set of all the edge objects present in the graph
        """
        edges = set()
        for adjacentVertex in self.adjacency_dictionary.values():
            edges.update(adjacentVertex.values())
        return edges

    def clear(self):
        """
        Clears all the parameters of the graph

        Performance: O(1)
        :return: None
        """
        self.adjacency_dictionary = {}
        self.size = 0
        self.paths = []
        self.positions = []

    def add_vertex(self, element):
        """
        Creates a new vertex object and adds it to the graph using the vertex_to_graph method.

        Performance: O(1)
        :param element: a string representing the element of the vertex
        :return: the vertex object which has been newly created
        """
        new_vertex = Vertex(element)
        self.vertex_to_graph(new_vertex)
        return new_vertex

    # TODO try for key error
    def vertex_to_graph(self, vertex):
        """
        Adds a vertex object to the graph by assigning a dictionary which will contain all adjacent vertices and edges.

        Performance: O(1)
        :param vertex: the vertex object that is to be added to the graph
        :return: None
        """
        self.adjacency_dictionary[vertex] = {}
        self.positions.append(vertex)
        self.size += 1

    # TODO try for key error
    def remove_vertex(self, vertex):
        """
        Delete the vertex object from the graph by removing it from the dictionaries of vertices which are adjacent.

        Performance: O(degree of vertex)
        :param vertex: the vertex object that is to be removed
        :return: None
        """
        for neighbour in self.neighbours(vertex):
            del self.adjacency_dictionary[neighbour][vertex]
        del self.adjacency_dictionary[vertex]
        self.positions.remove(vertex)

    def add_edge(self, first_vertex, second_vertex, element=None):
        """
        Create an edge object and add it to the graph using the edge_to_graph method.

        Performance: O(1)
        :param first_vertex: a vertex object that is to be an endpoint of the edge
        :param second_vertex: a vertex object that is to be an endpoint of the edge
        :param element: the element of the edge object
        :return: the edge object which has been newly added to the graph
        """
        new_edge = Edge(first_vertex, second_vertex, element)
        self.edge_to_graph(first_vertex, second_vertex, new_edge)
        return new_edge

    def edge_to_graph(self, first_vertex, second_vertex, edge):
        """
        Adds an edge object to the graph by adding it to the adjacency dictionary entries for both of its endpoints.

        Performance: O(1)
        :param first_vertex: one of the endpoints of the edge that is being added
        :param second_vertex: one of the endpoints of the edge that is being added
        :param edge: the edge that is to be added to the graph
        :return: None
        """
        # Try for key error
        self.adjacency_dictionary[first_vertex][second_vertex] = edge
        self.adjacency_dictionary[second_vertex][first_vertex] = edge

    def remove_edge(self, first_vertex, second_vertex):
        """
        Remove the edge that is found between the two given vertices.

        Performance: O(1)
        :param first_vertex:
        :param second_vertex:
        :return:
        """
        del self.adjacency_dictionary[first_vertex][second_vertex]
        del self.adjacency_dictionary[second_vertex][first_vertex]

    def neighbours(self, vertex):
        """
        Returns the vertices adjacent to the given vertex in the graph

        :param vertex: the vertex whose neighbours are to be found
        :return: a list of vertex objects adjacent to the vertex
        """
        return self.adjacency_dictionary[vertex].keys()

    def connecting_edges(self, vertex):
        """
        Returns the edges which are attached to the given vertex in the graph

        :param vertex: the vertex whose attached edges are to be found
        :return: a list of edge objects which are attached to the given vertex
        """
        return self.adjacency_dictionary[vertex].values()

    def degree(self, vertex):
        """
        Returns the degree of the given vertex which is the number of edges that are attached to the vertex

        :param vertex: the vertex object whose degree is to be returned
        :return: an integer which is the number of edges attached to the vertex
        """
        return len(self.adjacency_dictionary[vertex])

    def position_of_vertex(self, vertex):
        """
        Returns the position of the given vertex

        :param vertex: the vertex object whose position in the positions list is to be returned
        :return: an integer which is the position of the vertex in the graph
        """
        return self.positions.index(vertex)

    def vertex_from_position(self, position):
        """
        Tests if there is a vertex at that position and return it if there is

        :param position: an integer which is to be tested
        :return: a vertex if it had has the given position
        """
        if len(self.positions) > position:
            return self.positions[position]
        else:
            return False

    def endpoints_position(self, edge):
        """
        Gives the positions of the endpoints of the edge
        :return: a tuple of the positions of the endpoints
        """
        endpoints = edge.endpoints
        return self.position_of_vertex(endpoints[0]), self.position_of_vertex(endpoints[1])

    def contains_edge(self, first_vertex, second_vertex):
        """
        Test if there is an edge joining two vertices in the graph

        :param first_vertex: a vertex that is to be tested for a connection
        :param second_vertex: a vertex that is to be tested for a connection
        :return: an edge object which joins the two given vertices
        """
        if second_vertex in self.adjacency_dictionary[first_vertex]:
            return self.adjacency_dictionary[first_vertex][second_vertex]
        else:
            return False

    def find_all_paths(self):
        """
        Finds all the possible paths in a graph using a depth first search

        The depth first search is carried out starting from each vertex of the graph so all path combinations are found
        Algorithm structure from Handbook of Graph Theory, Gross & Yellen
        :return: a dictionary containing the strings representing the paths and their lengths as the value
        """
        completed = []      # The nodes which have acted as a root for the search
        all_paths = {}      # The dictionary of all the paths and their lengths (dict removes duplicate paths)
        for v in self.vertices():
            if v not in completed:
                for w in self.vertices():
                    w.visited = False       # Start search anew for each root
                path_stack = []             # Used to create the string of the path
                vertex_stack = []         # Used to create the list of vertices
                self._path_visit(v, path_stack, vertex_stack, all_paths)
                completed.append(v)
        return all_paths

    def _path_visit(self, v, path_stack, vertex_stack, all_paths):
        """
        The recursive element of the find_all_paths method.
        As each new path is found the string representing it is created using the path_stack.
        The vertices which are present in the path are stored along with the path as a tuple in self.paths
        :param v: the vertex object which is to be added to the path
        :param path_stack: The stack which is used to construct the path string
        :param vertex_stack: The stack which is used to store the vertices encountered in order
        :param all_paths: A dictionary which is used to store a non-duplicated list of paths along with their lengths
        :return: None
        """
        v.visited = True
        path_stack.append(v.element + '-')
        vertex_stack.append(v)
        for w in self.adjacency_dictionary[v]:
            if not w.visited:
                self._path_visit(w, path_stack, vertex_stack, all_paths)
        letters = ''.join(path_stack)
        path = letters[:-1]    # Remove final dash to make string more readable
        positions = list(vertex_stack)
        all_paths[path] = len(letters)/2
        self.paths.append((path, positions))
        # Once the vertex has been processed it is popped from the stack to create the string and to store the vertices
        path_stack.pop()
        vertex_stack.pop()

    def adjacency_dictionary_display(self):
        """
        Creates a string format of the adjacency dictionary of the given structure

        :return: string representation of the adjacency dictionary
        """
        string_list = []
        for vertex in self.adjacency_dictionary:
            string_list.append(str(vertex) + ':\n')
            for neighbour in self.neighbours(vertex):
                string_list.append(
                    '        ' + str(neighbour) + ': ' + str(self.adjacency_dictionary[vertex][neighbour]) + '\n')
        display_string = ''.join(string_list)
        return display_string
