import networkx as nx

class Graph:
    """ Take in a matrix and create a graph from it. The graph is initialized with a given number of rows and columns.

    Args:
    input_matrix (np.ndarray): The matrix to create the graph from.
    
    """
    def __init__(self, input_matrix, x_labels, y_labels):
        self.input_matrix = input_matrix
        self.x_labels = x_labels
        self.y_labels = y_labels
        
    
    def initialize_graph1(self):
        """ Initialize the graph with nodes and edges. """
        self.G = nx.Graph()
        self.G.add_nodes_from(self.x_labels)
        self.G.add_nodes_from(self.y_labels)
        #self.G.add_edges_from(np.argwhere(self.input_matrix == 1))
        # Add edges between nodes based on the input matrix
        for i in range(len(self.input_matrix)):
            for j in range(len(self.input_matrix[0])):
                if self.input_matrix[i][j] == 1:
                    self.G.add_edge(self.x_labels[i], self.y_labels[j])


    def draw_graph1(self, include_label:bool):
        """ Draw the graph using NetworkX. """
        self.initialize_graph1()
        # Draw the graph
        # blue for nodes in rows, red for nodes in columns
        node_color = ['blue' if node in self.x_labels else 'red' for node in self.G.nodes()] 
        # show nodes without any connections
        nx.draw(self.G, with_labels= True if include_label is True else False, node_color=node_color)
