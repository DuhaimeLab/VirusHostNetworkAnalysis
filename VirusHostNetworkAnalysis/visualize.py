import networkx as nx
import matplotlib.pyplot as plt
from sympy import deg

class Graph:
    """ Take in a matrix and create a graph from it. The graph is initialized with a given number of rows and columns.

    Args:
    input_matrix (np.ndarray): The matrix to create the graph from.
    
    """
    def __init__(self, input_matrix, x_labels, y_labels):
        self.input_matrix = input_matrix
        self.x_labels = x_labels
        self.y_labels = y_labels
        
        #print(type(self.input_matrix))
    
    def initialize_graph(self):
        """ Initialize the graph with nodes and edges. """
        self.G = nx.Graph()
        self.G.add_nodes_from(self.x_labels)
        self.G.add_nodes_from(self.y_labels)
        # self.G.add_edges_from(np.argwhere(self.input_matrix == 1))
        # Add edges between nodes based on the input matrix
        for i in range(len(self.input_matrix)):
            for j in range(len(self.input_matrix[1])):
                if self.input_matrix[i][j] == 1:
                    self.G.add_edge(self.x_labels[i], self.y_labels[j])


    def draw_graph(self, include_label:bool):
        """ Draw the graph using NetworkX. """
        plt.figure(figsize=(40, 30))
        self.initialize_graph()
        # Draw the graph. Blue for nodes in rows, red for nodes in columns
        node_color = ['blue' if node in self.x_labels else 'red' for node in self.G.nodes()]

        # Add one to that nodes with 0 degrees are visible
        node_size = [(self.G.degree(node)+1) * 50 for node in self.G.nodes()]

        # Set node size proportional to the degree of the node
        pos = nx.random_layout(self.G, seed=42)
        nx.draw(self.G, pos, with_labels= True if include_label is True else False, node_color=node_color,
                node_size = 100) 
        
        # Use gephi adjust graph for better visualization
    
    # Calculate the centrality of the graph
    def calculate_centrality(self):
        """ Calculate the centrality of the graph. """
        degree = nx.degree_centrality(self.G)
        degree_avg = sum(degree.values()) / len(degree)
        return degree_avg
    
    def degree_distribution(self, degree_seq):
        plt.figure(figsize=(10, 6))
        virus_degree = degree_seq[0]
        host_degree = degree_seq[1]

        # Plot histograms for virus and host degree distributions side by side)
        plt.hist(virus_degree, color='g', density=True)
        plt.title('Virus Degree Distribution')
        plt.xlabel('Degree')
        plt.ylabel('Frequency')
        plt.grid()
        plt.show()
    
        plt.hist(host_degree, color='b', density=True)
        plt.title('Host Degree Distribution')
        plt.xlabel('Degree')
        plt.ylabel('Frequency')
        plt.grid()
        plt.show()
       