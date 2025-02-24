import numpy as np
import networkx as nx

class ER:
    """ Class to create a random graph using the Erdős-Rényi model. 
    The graph is initialized with a given number of rows and columns, and a probability p for edge creation.
    
    Args:
    rows (int): Number of rows in the graph.
    columns (int): Number of columns in the graph.
    p (float): Probability of edge creation between nodes.

    """
    def __init__(self, rows:int, columns:int, p:float):
        self.rows = list(range(0, rows))
        self.columns = list(range(rows, rows + columns))
        self.p = p
        self.matrix_rand = np.zeros((len(self.rows), len(self.columns)), dtype=bool)

    def initialize_graph(self):
        """ Initialize the graph with nodes and edges. """
        self.G = nx.Graph()
        self.G.add_nodes_from(self.rows)
        self.G.add_nodes_from(self.columns)

    def fill_ER_graph(self):
        """ Create a graph with n nodes and random edges between them."""
       # self.initialize_graph()
        # Iterate through all pairs of nodes in the graph
        for row in self.rows:
            for col in self.columns:
                p_rand = np.random.rand()
                # if the random probability is less than p, add an edge between the nodes
                if p_rand < self.p:
                    #self.G.add_edge(row, col)
                    self.matrix_rand[row][col - len(self.rows)] = 1
        print(self.matrix_rand)
        return self.matrix_rand

    # def calculate_centrality(self):
    #     self.betweenness = nx.betweenness_centrality(self.G)
    #     self.closeness = nx.closeness_centrality(self.G)
    #     self.degree = nx.degree_centrality(self.G)
    #     self.eigenvector = nx.eigenvector_centrality(self.G)
    #     self.pagerank = nx.pagerank(self.G)
    #     print(self.degree)


# Configuration Model 
#class CM:

