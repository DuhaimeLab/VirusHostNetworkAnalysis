from networkx import eigenvector_centrality, random_clustered_graph
import numpy as np
import networkx as nx
import igraph as ig
import random
import matplotlib.pyplot as plt

from sympy import degree

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
        return self.matrix_rand

    def calculate_centrality(self):
        # self.betweenness = nx.betweenness_centrality(self.G)
        # self.closeness = nx.closeness_centrality(self.G)
        self.degree = nx.degree_centrality(self.G)
        # self.eigenvector = nx.eigenvector_centrality(self.G)
        # self.pagerank = nx.pagerank(self.G)
        #print(self.degree)
        # find average centrality
        avg_degree_centrality = sum(self.degree.values()) / len(self.degree)
        print(f"Average Degree Centrality: {avg_degree_centrality:.4f}")
        print("done")


# Configuration Model 
class CM:
    """ Class to create a random graph using the Configuration Model. 
    The graph is initialized with a given number of rows and columns, and a probability p for edge creation.
    
    Args:
    rows (int): Number of rows in the graph.
    columns (int): Number of columns in the graph.
    p (float): Probability of edge creation between nodes.

    """
    def __init__(self, matrix_vhip, edge_list:list):
        """ Initialize the Configuration Model with a given matrix and edge list."""
        # self.rows = list(range(0, rows))
        # self.columns = list(range(rows, rows + columns))
        # self.p = p
        self.matrix_vhip = matrix_vhip
        self.edge_list = edge_list

    def initialize_graph(self):
        """ Initialize the graph using the actual data."""
        G = nx.Graph()
        G.add_edges_from(self.edge_list)

    
    def find_candidates(self):
        """ Randomly selects two viruses and finds candidate edges to swap between them."""
        failed_runs = 0
        successful_runs = 0
        # number of viruses/columns
        num_viruses = len(self.matrix_vhip)

        # pick random virus indices
        self.index1 = random.randint(0, num_viruses - 1)  # random index for virus 1
        self.index2 = random.randint(0, num_viruses - 1)  # random index for virus 2
        while self.index1 == self.index2:  # ensure they are different
            self.index2 = random.randint(0, num_viruses - 1)
        print(f"Selected virus indices: {self.index1}, {self.index2}")

        # store the rows corresponding to the random virus indices
        virus_array1 = self.matrix_vhip[self.index1]  # get the row for random virus index 1
        virus_array2 = self.matrix_vhip[self.index2]  # get the row for random virus index 2

        # initialize empty lists for candidate edges
        virus1_candidates = []
        virus2_candidates = []

        # loop
        for counter, (v1, v2) in enumerate(zip(virus_array1, virus_array2)):
            if v1 == v2:
                continue
            else:
                if v1 == 0:
                    virus2_candidates.append(counter)
                else: 
                    virus1_candidates.append(counter)
        
        # Check if empty
        if not virus1_candidates or not virus2_candidates:
            failed_runs += 1
            print("failed to find candidates: one of the lists is empty.")
        else:
            print(f"Virus 1 candidates: {virus1_candidates}")
            print(f"Virus 2 candidates: {virus2_candidates}")
            successful_runs += 1
            return virus1_candidates, virus2_candidates
        
    def update_matrix(self, virus1_candidates:list, virus2_candidates:list):
        """ Update the matrix with the new edges randomly swapped in the previous step. """
        # randomly pick one candidate from each list
        random_col1 = random.choice(virus2_candidates)
        random_col2 = random.choice(virus1_candidates)

        # Update values
        self.matrix_vhip[self.index1][random_col1] = 1
        self.matrix_vhip[self.index1][random_col2] = 0
        self.matrix_vhip[self.index2][random_col1] = 0
        self.matrix_vhip[self.index2][random_col2] = 1

    def run_config_model(self):
        # Check that function ran successfully
        vals = self.find_candidates()
        if vals is not None:
            virus1_candidates, virus2_candidates = vals
            self.update_matrix(virus1_candidates, virus2_candidates)
            print(self.matrix_vhip)

    def bootstrap_stats(self, iterations):
        """ Run the Configuration Model for a number of iterations to generate random graphs. """
        for i in range(iterations):
            self.run_config_model()
        
        # create edge list from the updated matrix
        edge_list = []
        for i in range(len(self.matrix_vhip)):
            for j in range(len(self.matrix_vhip[0])):
                if self.matrix_vhip[i][j] == 1:
                    edge_list.append((i + len(self.matrix_vhip[0]), j))
        print(f"Generated edge list after {iterations} iterations: {edge_list}")

        G = nx.Graph()
        G.add_edges_from(edge_list)
        # plot degree distribution
        self.degree_sequence = sorted([d for n, d in G.degree()], reverse=True)


    
    def iterations(self, iterations):
        stats = []
        for i in range(iterations):
            avg_betweenness = self.bootstrap_stats(iterations)
            stats.append(avg_betweenness)

        # plot stats in a line graph where x-axis is the iteration number and y-axis is the average betweenness centrality
        plt.figure(figsize=(10, 6))
        plt.hist(self.degree_sequence, bins= 30, color='b')
        plt.title('Degree Distribution')
        plt.xlabel('Degree')
        plt.ylabel('Frequency')
        plt.grid()
        plt.show()
        return stats

        



        

    


    

