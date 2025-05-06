import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import matplotlib.colors as mcol
import networkx as nx
from multiprocessing import Pool
from typing import List, Tuple
from tqdm import tqdm
from networkx.algorithms import community


class PredictionMatrix:
    """ Class to make a matrix from the input file. Can make square or rectangle matrices. Can also be used to make probabilty matrix.
    
    Args:
    vhip_predictions_file (str): Path to the input file containing virus-host interaction predictions.
    
    """

    def __init__(self, vhip_predictions_file:str, probability:bool = False):
        #self.file = vhip_predictions_file
        #check if this file path exists
        if not os.path.exists(vhip_predictions_file):
            raise FileNotFoundError(f"The file {vhip_predictions_file} does not exist.")
        else:
            self.virus_host = pd.read_csv(vhip_predictions_file, sep='\t')
            self.error_check()
            self.predictions = self.virus_host[self.virus_host['Predictions'] == 1]
            self.title = vhip_predictions_file.split('/')[-1].split('.')[0]
            self.probability = probability

    def error_check(self):
        """ Check if the input file is in the correct format. """
        # Check if the file has the three columns needed for visualization
        # The files need pairs, Predicitons, and either InfProbabilities OR Scores
        # the file does not need to have all three columns, but it needs to have at least one of the last two
        # Check if the file has the required columns
        required_columns = ['pairs', 'Predictions']
        if not all(col in self.virus_host.columns for col in required_columns):
            raise ValueError(f"The input file must contain the following columns: {', '.join(required_columns)}")
        # Check if the file has either InfProbabilities or Scores
        if not any(col in self.virus_host.columns for col in ['InfProbabilities', 'Scores']):
            raise ValueError("The input file must contain either the InfProbabilities or Scores column.")
        # Check if the pairs column is in the correct format
        if not all(self.virus_host['pairs'].str.contains(':')):
            raise ValueError("The pairs column is not in the correct format. It should be in the format 'virus:host'")

    def get_unique_virus_host(self):
        """Find the unique virus names and host names in the dataset. The first column contains a virus name and host name, which are separated by a colon."""
        # The unique viruses will make of the rows.
        self.rows = self.virus_host['pairs'].str.split(':').str[0].unique()
        # The unique hosts will make up the columns.
        self.columns = self.virus_host['pairs'].str.split(':').str[1].unique()

    def initialize_matrix(self):
        """"Create a matrix full of zeros and make a list of rows and column labels.
        Args:
        matrix_type (str): Type of matrix to create, which determines if the values are boolean or float.
        
        """""
        # self.rows = self.unique_viruses
        # self.columns = self.unique_hosts
        self.virus_host_array = np.zeros((len(self.rows), len(self.columns)), dtype= bool if self.probability is False else float)
    
    def fill_matrix(self):
        """ Fill the prediction matrix with 1s and 0s or fill the probability matrix with the InfProbabilities from the dataset."""
        if self.probability is False:
            # Only loops through the subset of data where predicitons == 1
            for _, row in self.predictions.iterrows():
                virus, host = row['pairs'].split(':')
                # Find the index of the virus (row) and host (col) in the matrix
                virus_index = np.where(self.rows == virus)[0][0]
                host_index = np.where(self.columns == host)[0][0]
                # Fill the matrix with 1 if matrix_type == 'predictions'
                self.virus_host_array[virus_index][host_index] = 1
        else:
            for _, row in self.virus_host.iterrows():
                virus, host = row['pairs'].split(':')
                # Find the index of the virus (row) and host (col) in the matrix
                virus_index = np.where(self.rows == virus)[0][0]
                host_index = np.where(self.columns == host)[0][0]
                # Fill the matrix with InfProbabilities or Scores depending on the column name
                if 'InfProbabilities' in row:
                    self.virus_host_array[virus_index][host_index] = row['InfProbabilities']
                elif 'Scores' in row:
                    self.virus_host_array[virus_index][host_index] = row['Scores']
                
                #print(row['InfProbabilities'])

    def make_rectangular_matrix(self):
        """ Call all functions to make a rectangular matrix"""
        self.get_unique_virus_host()
        self.initialize_matrix()
        self.fill_matrix()
        return self.virus_host_array

    def save_matrix(self):
        """ Save the matrix to a csv file. """
        np.savetxt('virus_host.csv', self.virus_host_array, delimiter=',')
        print('virus_host.csv saved')

    # def create_edge_list(self):
    #     """ Create an edge list from the predictions column of the dataset."""""
    #     edge_list = []
    #     for _, row in self.predictions.iterrows():
    #         virus, host = row['pairs'].split(':')
    #         edge_list.append((str(virus), str(host)))
    #     print("Edge list created with {} edges.".format(len(edge_list)))
    #     # save edge list to a csv file
    #     edge_df = pd.DataFrame(edge_list, columns=['Virus', 'Host'])
    #     edge_df.to_csv('edge_list.csv', index=False)
    #     return edge_list



# class Calculations:
#     """ Class to calculate nestedness, degree, modularity, etc (to be expanded later).
    
#     Args:
#     mat (2d array): Matrix to be used for calculations.
#     sorted (bool): If the matrix is already sorted or not. Default is False.
    
#     """

#     def __init__(self, mat, sorted:bool = False):
#         # Initialize the class with a sorted adjacency matrix
#         self.mat = np.array(mat)
#         self.sorted = sorted
#         if sorted is False:
#             self.sort()

#     def sort_rows_cols(self, axis:int):
#         """ Find the sum of each row and move rows with the highest sum to the top of the matrix. """
#         counts = np.sum(self.mat, axis=axis)
#         # Sort inidices based on the counts
#         sorted_indices = np.argsort(counts)[::-1]
#         # Sort the matrix and based on the sorted indices
#         if axis == 1:
#             self.mat = self.mat[sorted_indices]
#         elif axis == 0:
#             self.mat= self.mat[:, sorted_indices]
#     def sort(self):
#         self.sort_rows_cols(1)
#         self.sort_rows_cols(0)


#     def nestedness_rows(self, pair):
#         """Calculate nestedness for rows using the NODF algorithm for the array."""
#         # generate list of rows to compare
#         pair1 = self.mat[pair[0],]
#         pair2 = self.mat[pair[1],]
#         N_row = self.compare(pair1, pair2)
#         #print(pair[0], pair[1], N_row)
#         #print(N_row)
#         return N_row

#     def nestedness_cols(self, pair):
#         """Calculate nestedness for cols using the NODF algorithm for the array."""
#         pair1 = self.mat[:, pair[0]]
#         pair2 = self.mat[:, pair[1]]
#         N_col = self.compare(pair1, pair2)  # pyright: ignore
#         return N_col

#     def pairs(self, axis: int = 0) -> list[tuple[int, int]]:
#         """Determine all possible i-j pairs.

#         Args:
#             axis (int): Axis to be used when determining all pairs.
#         """
#         lst: List[tuple[int, int]] = []
#         for i in range(0, self.mat.shape[axis]):
#             for j in range(i + 1, self.mat.shape[axis]):
#                 lst.append((i, j))
#         return lst
#         # equation for number of pairs


#     def compare(self, x: list[int], y: list[int]) -> float:
#         """Compare two lists containing 0 and 1.

#         Args:
#             x (list[int]): first list
#             y (list[int]): second list
#         """
#         if sum(x) <= sum(y):
#             val = 0
#         elif sum(y) == 0:
#             val = 0
#         else:
#             counter = 0
#             total = 0
#             for i, j in zip(x, y):
#                 if i == 1 and j == 1:
#                     counter += 1
#                     total += 1
#                 elif i == 0 and j == 1:
#                     total += 1
#             val = counter / total
#         return val * 100

#     def run_parallel(self, num_procs: int = 6):
#         """Run multiple process of the compute_feature method.

#         This significantly improves run time by using multiple CPU cores.

#         Args:
#             num_procs (int): Number of core to be used.
#         """

#         with Pool(num_procs) as pool:
#             #nrow = pool.map(self.nestedness_rows, self.pairs(axis=0))
#             with tqdm(total=len(self.pairs(axis=0)) + len(self.pairs(axis=1)), desc="Calculating nestedness", colour="green") as pbar:
#                 nrow = []
#                 for result in pool.imap_unordered(self.nestedness_rows, self.pairs(axis=0)):
#                     nrow.append(result)
#                     pbar.update(1)

#                 ncol = []
#                 for result in pool.imap_unordered(self.nestedness_cols, self.pairs(axis=1)):
#                     ncol.append(result)
#                     pbar.update(1)

#         nrow = sum(nrow) / (len(self.mat) * (len(self.mat) - 1) / 2)
#         ncol = sum(ncol) / (len(self.mat[0]) * (len(self.mat[0]) - 1) / 2)
#         print(nrow, ncol)
#         nodf = (nrow + ncol) / 2
#         return nodf
    
#     def calculate_degree(self):
#         """ Calculate the degree of each node in the graph. """
#         host_degrees = []
#         for col in range(len(self.mat[0])):
#             total = 0
#             for row in range(len(self.mat)):
#                 if self.mat[row][col] == 1:
#                     total += 1
#             host_degrees.append(total)
        
#         virus_degrees = []
#         for row in range(len(self.mat)):
#             total = 0
#             for col in range(len(self.mat[0])):
#                 if self.mat[row][col] == 1:
#                     total += 1
#             virus_degrees.append(total)

#         return [virus_degrees, host_degrees]
    

#     # Calculate modularity
#     def calculate_modularity(self):
#         """ Calculate the modularity of the graph. """
#         G = nx.from_numpy_array(self.mat)
#         partition = community.greedy_modularity_communities(G)
#         modularity = community.modularity(G, partition)
#         print("Modularity: ", modularity)
#         return modularity