import numpy as np
import networkx as nx
import random
import matplotlib.pyplot as plt

class ER:
    """ Class to create a random graph using the Erdős-Rényi model. 
    The graph is initialized with the PredictionMatrix, which has a given number of rows and columns, and a probability p for edge creation.
    
    Args:
    PredictionMatrix (PredictionMatrix): A class that contains the true virus-host interaction matrix, row names, and column names.
    p (float): Probability of edge creation between nodes.
    """
    def __init__(self, PredictionMatrix, p:float):
        self.rows = PredictionMatrix.rows
        self.columns = PredictionMatrix.columns
        self.p = p
        # Initialize the ER array to be the same size as the prediction matrix, fill with all zeros
        self.virus_host_array = np.zeros((len(self.rows), len(self.columns)), dtype=bool)

    def fill_ER_graph(self):
        """ Create a graph with n nodes and random edges between them."""
        # Iterate through all pairs of nodes in the graph
        for row in range(0, len(self.rows)):
            for col in range(0, len(self.columns)):
                p_rand = np.random.rand()
                # if the random probability is less than p, add an edge between the nodes
                if p_rand < self.p:
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


class ConfigurationModel:
    """ Class to create a random graph using the Configuration Model. 
     The graph is initialized with the PredictionMatrix, which has a given number of rows and columns,
    
    Args:
    PredictionMatrix (PredictionMatrix): A class that contains the true virus-host interaction matrix, row names, and column names.
    """
    def __init__(self, PredictionMatrix):
        """ Initialize the Configuration Model using the given prediction matrix. """
        self.virus_host_array= PredictionMatrix.virus_host_array
        self.rows = PredictionMatrix.rows
        self.columns = PredictionMatrix.columns
  
    def find_candidates(self):
        """ Randomly selects two viruses and finds candidate edges to swap between them."""
        # find number of viruses/columns
        num_viruses = len(self.virus_host_array)

        # pick random virus indices
        self.index1 = random.randint(0, num_viruses - 1)  # random index for virus 1
        self.index2 = random.randint(0, num_viruses - 1)  # random index for virus 2
        while self.index1 == self.index2:  # ensure they are different
            self.index2 = random.randint(0, num_viruses - 1)

        # store the rows corresponding to the random virus indices
        virus_array1 = self.virus_host_array[self.index1]  # get the row for random virus index 1
        virus_array2 = self.virus_host_array[self.index2]  # get the row for random virus index 2

        # initialize empty dictionaries for candidate edges
        virus1_candidates = {}
        virus2_candidates = {}

        # find candidate edges to swap
        # for counter, (v1, v2) in enumerate(zip(virus_array1, virus_array2)):
        #     if v1 == v2:
        #         continue
        #     else:
        #         if v1 == 0:
        #             virus2_candidates[counter] = self.index1
        #         else: 
        #             virus1_candidates[counter] = self.index2

        for counter, (v1, v2) in enumerate(zip(virus_array1, virus_array2)):
            if v1 != v2:
                if v1 == 0:
                    virus2_candidates[counter] = self.index1
                else:
                    virus1_candidates[counter] = self.index2

        # Check if either candidate dictionary is empty
        # If one or both dictionaries are empty, the swap cannot be performed
        if not bool(virus1_candidates) or not bool(virus2_candidates):
            self.failed_runs += 1
        else:
            self.successful_runs += 1
            return virus1_candidates, virus2_candidates
        

    def update_matrix(self, virus1_candidates, virus2_candidates):
        """ Update the matrix with the new edges randomly swapped in the previous step.
         
        Args:
        virus1_candidates (dict): Dictionary of candidates for virus 1. Key is the column index, value is the row index.
        virus2_candidates (dict): Dictionary of candidates for virus 2. Key is the column index, value is the row index.
        """
        # randomly pick one candidate from each list
        random_col1 = random.choice(list(virus2_candidates.keys()))
        random_col2 = random.choice(list(virus1_candidates.keys()))

        # Update values in the matrix
        self.virus_host_array[self.index1][random_col1], self.virus_host_array[self.index1][random_col2] = 1, 0
        self.virus_host_array[self.index2][random_col1], self.virus_host_array[self.index2][random_col2] = 0, 1

    def run_config_model(self):
        """ Run all steps of the Configuration Model for one swap. """
        # Check that function ran successfully
        candidates = self.find_candidates()
        if candidates:
            # If candidates are found, update the matrix
            self.update_matrix(candidates[0], candidates[1])

    def bootstrap_stats(self, swaps):
        """ Run the Configuration Model for a number of iterations to generate a random matrix. 
        Args:
        swaps (int): Number of successful swaps to make.
        """
        self.failed_runs = 0
        self.successful_runs = 0

        # run the configuration model for a number of successful swaps
        while self.successful_runs < swaps:
            self.run_config_model()

class RandomShuffle:
    """ Class to create a random graph. 
    The graph is initialized with a given number of rows and columns.
    
    Args:
    PredictionMatrix (PredictionMatrix): A class that contains the true virus-host interaction matrix, row names, and column names.
    """
    def __init__(self, PredictionMatrix):
        self.virus_host_array= PredictionMatrix.virus_host_array
        self.rows = PredictionMatrix.rows
        self.columns = PredictionMatrix.columns

    def shuffle_ones(self):
        """ Shuffle the 1s in the matrix. """
        # Get sum of the matrix
        sum_matrix = np.sum(self.virus_host_array)
        # Randomly select sum_matrix number of spots in the matrix
        random_spots = random.sample(range(len(self.virus_host_array.flatten())), sum_matrix)
        # Create a new matrix with all 0s
        new_matrix = np.zeros_like(self.virus_host_array)
        # Set the selected spots to 1
        for spot in random_spots:
            row = spot // self.virus_host_array.shape[1]
            col = spot % self.virus_host_array.shape[1]
            new_matrix[row][col] = 1
        # Return the new matrix
        print(self.virus_host_array.shape)
        self.virus_host_array = new_matrix
        return new_matrix
    
    

        



        

    


    

