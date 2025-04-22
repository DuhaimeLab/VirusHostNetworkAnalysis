from networkx import eigenvector_centrality
import numpy as np
import networkx as nx
import random
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
import seaborn as sns


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

    # Code shows that even when sorted, random networks will be moslty evenly spread
    # def sort_rows_cols(self, axis:int):
    #     """ Find the sum of each row and move rows with the highest sum to the top of the matrix. """
    #     counts = np.sum(self.matrix_rand, axis=axis)
    #     # Sort inidices based on the counts
    #     sorted_indices = np.argsort(counts)[::-1]
    #     # Sort the matrix and row names based on the sorted indices
    #     if axis == 1:
    #         self.rows = np.array(self.rows)[sorted_indices].tolist()
    #         self.matrix_rand = self.matrix_rand[sorted_indices]
    #     elif axis == 0:
    #         self.columns = np.array(self.columns)[sorted_indices].tolist()
    #         self.matrix_rand = self.matrix_rand[:, sorted_indices]
    # def sort_matrix(self):
    #     self.sort_rows_cols(1)
    #     self.sort_rows_cols(0)


    def fill_ER_graph(self):
        """ Create a graph with n nodes and random edges between them."""
        # Iterate through all pairs of nodes in the graph
        for row in range(0, len(self.rows)):
            for col in range(0, len(self.columns)):
                p_rand = np.random.rand()
                # if the random probability is less than p, add an edge between the nodes
                if p_rand < self.p:
                    #self.G.add_edge(row, col)
                    self.matrix_rand[row][col] = 1
        return self.matrix_rand

    def create_edge_list(self):
        """ Create an edge list from the matrix. """
        edge_list = []
        for i in range(len(self.matrix_rand)):
            for j in range(len(self.matrix_rand[0])):
                if self.matrix_rand[i][j] == 1:
                    edge_list.append((i + len(self.matrix_rand[0]), j))
        self.G = nx.Graph()
        self.G.add_edges_from(edge_list)
        #print(f"Generated edge list: {edge_list}")
        return edge_list
    
    def draw_graph(self, include_label:bool):
        """ Draw the graph using NetworkX. """
        plt.figure(figsize=(40, 30))
        # Draw the graph. Blue for nodes in rows, red for nodes in columns
        node_color = ['blue' if node in self.rows else 'red' for node in self.G.nodes()]

        # # Add one to that nodes with 0 degrees are visible
        # node_size = [(self.G.degree(node)+1) * 50 for node in self.G.nodes()]

        # Set node size proportional to the degree of the node
        pos = nx.random_layout(self.G, seed=42)
        nx.draw(self.G, pos, with_labels= True if include_label is True else False, node_color=node_color,
                node_size = 100) 
        
    def calculate_degree(self):
         # degree sequence for just the hosts
        host_degrees = []
        for col in range(len(self.matrix_rand[0])):
            total = 0
            for row in range(len(self.matrix_rand)):
                if self.matrix_rand[row][col] == 1:
                    total += 1
            host_degrees.append(total)
        
        virus_degrees = []
        for row in range(len(self.matrix_rand)):
            total = 0
            for col in range(len(self.matrix_rand[0])):
                if self.matrix_rand[row][col] == 1:
                    total += 1
            virus_degrees.append(total)

        return [virus_degrees, host_degrees]

    def plot_heatmap(self, matrix_type:str):
        cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",["r", "white", "b"])
        sns.heatmap(self.matrix_rand, cmap= mcol.LinearSegmentedColormap.from_list("MyCmapName",["white", "cadetblue"]) if matrix_type == 'prediction' else cm1)
        plt.gcf().set_size_inches(7, 14)
        plt.xlabel("Hosts")
        plt.ylabel("Viruses")
        # Replace 'your_file_name_here' with the actual file name or pass it as an argument to the class
        matrix_title = 'your_file_name_here'.replace('Sample_Input/', '')
        plt.title(matrix_title.split('_')[0] + ' ' + matrix_type)
        plt.savefig('Heatmaps/Heatmap_' + matrix_title.split('_')[0] + '_' +  matrix_type + '.png')
        plt.show()



# Configuration Model 
class CM:
    """ Class to create a random graph using the Configuration Model. 
    The graph is initialized with a given number of rows and columns, and a probability p for edge creation.
    
    Args:
    rows (int): Number of rows in the graph.
    columns (int): Number of columns in the graph.
    p (float): Probability of edge creation between nodes.

    """
    def __init__(self, matrix_vhip):
        """ Initialize the Configuration Model with a given matrix and edge list."""
        # self.rows = list(range(0, rows))
        # self.columns = list(range(rows, rows + columns))
        # self.p = p
        self.matrix_vhip = matrix_vhip
  
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
        #print(f"Selected virus indices: {self.index1}, {self.index2}")

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
            #print("failed to find candidates: one of the lists is empty.")
        else:
            #print(f"Virus 1 candidates: {virus1_candidates}")
            #print(f"Virus 2 candidates: {virus2_candidates}")
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
            #print(self.matrix_vhip)

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

        # degree sequence for just the hosts
        host_degrees = []
        for col in range(len(self.matrix_vhip[0])):
            total = 0
            for row in range(len(self.matrix_vhip)):
                if self.matrix_vhip[row][col] == 1:
                    total += 1
            host_degrees.append(total)
        
        virus_degrees = []
        for row in range(len(self.matrix_vhip)):
            total = 0
            for col in range(len(self.matrix_vhip[0])):
                if self.matrix_vhip[row][col] == 1:
                    total += 1
            virus_degrees.append(total)

        return [virus_degrees, host_degrees]
    
    # def sort_rows_cols(self, axis:int):
    #     """ Find the sum of each row and move rows with the highest sum to the top of the matrix. """
    #     counts = np.sum(self.matrix_vhip, axis=axis)
    #     # Sort inidices based on the counts
    #     sorted_indices = np.argsort(counts)[::-1]
    #     # Sort the matrix and based on the sorted indices
    #     if axis == 1:
    #         self.matrix_vhip = self.matrix_vhip[sorted_indices]
    #     elif axis == 0:
    #         self.matrix_vhip = self.matrix_vhip[:, sorted_indices]

    # def sort_matrix(self):
    #     self.sort_rows_cols(1)
    #     self.sort_rows_cols(0)

    def iterations(self, bootstraps, iterations):
        for i in range(iterations):
           [virus_degrees, host_degrees] = self.bootstrap_stats(bootstraps)
           #self.matrix_vhip
        return [virus_degrees, host_degrees]
    

    

        



        

    


    

