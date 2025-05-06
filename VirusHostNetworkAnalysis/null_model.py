import numpy as np
import networkx as nx
import random
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
from pytest import skip
import seaborn as sns
from streamlit import success
from tqdm import tqdm
from multiprocessing import Pool

class ER:
    """ Class to create a random graph using the Erdős-Rényi model. 
    The graph is initialized with a given number of rows and columns, and a probability p for edge creation.
    
    Args:
    PredictionMatrix (PredictionMatrix): A class that contains the true virus-host interaction matrix, row names, and column names.
    p (float): Probability of edge creation between nodes.

    """
    def __init__(self, PredictionMatrix, p:float):
        self.rows = PredictionMatrix.rows
        self.columns = PredictionMatrix.columns
        self.p = p
        self.virus_host_array = np.zeros((len(self.rows), len(self.columns)), dtype=bool)

    def fill_ER_graph(self):
        """ Create a graph with n nodes and random edges between them."""
        # Iterate through all pairs of nodes in the graph
        for row in range(0, len(self.rows)):
            for col in range(0, len(self.columns)):
                p_rand = np.random.rand()
                # if the random probability is less than p, add an edge between the nodes
                if p_rand < self.p:
                    #self.G.add_edge(row, col)
                    self.virus_host_array[row][col] = 1
        return self.virus_host_array

    def create_edge_list(self):
        """ Create an edge list from the matrix. """
        edge_list = []
        for i in range(len(self.virus_host_array)):
            for j in range(len(self.virus_host_array[0])):
                if self.virus_host_array[i][j] == 1:
                    edge_list.append((i + len(self.virus_host_array[0]), j))
        self.G = nx.Graph()
        self.G.add_edges_from(edge_list)
        #print(f"Generated edge list: {edge_list}")
        return edge_list
    
    def draw_graph(self, include_label:bool):
        """ Draw the graph using NetworkX. """
        plt.figure(figsize=(40, 30))
        # Draw the graph. Blue for nodes in rows, red for nodes in columns
        node_color = ['blue' if node in self.rows else 'red' for node in self.G.nodes()]

        # Set node size proportional to the degree of the node
        pos = nx.random_layout(self.G, seed=42)
        nx.draw(self.G, pos, with_labels= True if include_label is True else False, node_color=node_color,
                node_size = 100) 


# Configuration Model 
class ConfigurationModel:
    """ Class to create a random graph using the Configuration Model. 
    The graph is initialized with a given number of rows and columns, and a probability p for edge creation.
    
    Args:
    PredictionMatrix (PredictionMatrix): A class that contains the true virus-host interaction matrix, row names, and column names.
    """
    def __init__(self, PredictionMatrix):
        """ Initialize the Configuration Model with a given matrix and edge list."""
        self.virus_host_array= PredictionMatrix.virus_host_array
        self.rows = PredictionMatrix.rows
        self.columns = PredictionMatrix.columns
  
    def find_candidates(self):
        """ Randomly selects two viruses and finds candidate edges to swap between them."""
        # number of viruses/columns
        num_viruses = len(self.virus_host_array)

        # pick random virus indices
        self.index1 = random.randint(0, num_viruses - 1)  # random index for virus 1
        self.index2 = random.randint(0, num_viruses - 1)  # random index for virus 2
        while self.index1 == self.index2:  # ensure they are different
            self.index2 = random.randint(0, num_viruses - 1)
        #print(f"Selected virus indices: {self.index1}, {self.index2}")

        # store the rows corresponding to the random virus indices
        virus_array1 = self.virus_host_array[self.index1]  # get the row for random virus index 1
        virus_array2 = self.virus_host_array[self.index2]  # get the row for random virus index 2

        # initialize empty lists for candidate edges
        virus1_candidates = {}
        virus2_candidates = {}

        for counter, (v1, v2) in enumerate(zip(virus_array1, virus_array2)):
            if v1 == v2:
                continue
            else:
                if v1 == 0:
                    virus2_candidates[counter] = self.index1
                else: 
                    virus1_candidates[counter] = self.index2

        # Check if empty
        if not bool(virus1_candidates) or not bool(virus2_candidates):
            self.failed_runs += 1
            #print("failed to find candidates: one of the lists is empty.")
        else:
            #print(f"Virus 1 candidates: {virus1_candidates}")
            #print(f"Virus 2 candidates: {virus2_candidates}")
            self.successful_runs += 1
            return virus1_candidates, virus2_candidates
        

    def update_matrix(self, virus1_candidates, virus2_candidates):
        """ Update the matrix with the new edges randomly swapped in the previous step. """
        # randomly pick one candidate from each list
        random_col1 = random.choice(list(virus2_candidates.keys()))
        random_col2 = random.choice(list(virus1_candidates.keys()))

        # Update values
        self.virus_host_array[self.index1][random_col1] = 1
        self.virus_host_array[self.index1][random_col2] = 0
        self.virus_host_array[self.index2][random_col1] = 0
        self.virus_host_array[self.index2][random_col2] = 1

    def run_config_model(self):
        # Check that function ran successfully
        candidates = self.find_candidates()
        if candidates:
            self.update_matrix(candidates[0], candidates[1])

    def bootstrap_stats(self, iterations):
        """ Run the Configuration Model for a number of iterations to generate random graphs. 
        Args:
        iterations (int): Number of iterations to run the model.
        """
        self.failed_runs = 0
        self.successful_runs = 0

        # run the configuration model for a number of iterations
        #with tqdm(total=iterations, desc="Swapping edges", colour="green") as pbar:
        while self.successful_runs < iterations:
            self.run_config_model()
        #pbar.update(1)
        #print("Successful runs: ", self.successful_runs, "Failed runs: ", self.failed_runs)
    

    

        



        

    


    

