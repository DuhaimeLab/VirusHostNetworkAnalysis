import numpy as np
import networkx as nx

## ER Class
class ER:
    def __init__(self, rows:int, columns:int, p:float):
        self.rows = list(range(1, rows))
        self.columns = list(range(rows, rows + columns))
        self.p = p

    def initialize_graph(self):
        self.G = nx.Graph()
        self.G.add_nodes_from(self.rows)
        self.G.add_nodes_from(self.columns)
        print(self.rows)
        print(self.columns)

    # Create a graph with n nodes and random edges between them
    def fill_ER_graph(self):
        # loop through the nodes, rows then columns
        for row in self.rows:
            for col in self.columns:
                if row != col:
                    # set p to a random value between 0 and 1
                    p_rand = np.random.rand()
                    # if p is less than 0.1, add an edge between the nodes
                    if p_rand < self.p:
                        self.G.add_edge(row, col)

    def draw_graph(self):
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