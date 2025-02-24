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
        self.rows = list(range(1, rows))
        self.columns = list(range(rows, rows + columns))
        self.p = p

    def initialize_graph(self):
        """ Initialize the graph with nodes and edges. """
        self.G = nx.Graph()
        self.G.add_nodes_from(self.rows)
        self.G.add_nodes_from(self.columns)

    def fill_ER_graph(self):
        """ Create a graph with n nodes and random edges between them."""

        # Iterate through all pairs of nodes in the graph
        for row in self.rows:
            for col in self.columns:
                if row != col:
                    # set p to a random value between 0 and 1
                    p_rand = np.random.rand()
                    # if the random probability is less than p, add an edge between the nodes
                    if p_rand < self.p:
                        self.G.add_edge(row, col)

    def draw_graph(self):
        """ Draw the graph using NetworkX. """
        self.initialize_graph()
        self.fill_ER_graph()
        # Draw the graph
        # blue for nodes in rows, red for nodes in columns
        node_color = ['blue' if node in self.rows else 'red' for node in self.G.nodes()] 
        # show nodes without any connections
        nx.draw(self.G, with_labels=True, node_color=node_color)

    def calculate_centrality(self):
        self.betweenness = nx.betweenness_centrality(self.G)
        self.closeness = nx.closeness_centrality(self.G)
        self.degree = nx.degree_centrality(self.G)
        self.eigenvector = nx.eigenvector_centrality(self.G)
        self.pagerank = nx.pagerank(self.G)
        print(self.degree)


# Configuration Model 
