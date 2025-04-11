import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import matplotlib.colors as mcol
import networkx as nx
import igraph as ig
import rustworkx as rx
from rustworkx.visualization import mpl_draw

class PredictionMatrix:
    """ Class to make a matrix from the input file. Can make square or rectangle matrices. Can also be used to make probabilty matrix.
    
    Args:
    vhip_predictions_file (str): Path to the input file containing virus-host interaction predictions.
    
    """

    def __init__(self, vhip_predictions_file:str):
        self.file = vhip_predictions_file
        #check if this file path exists
        if not os.path.exists(self.file):
            raise FileNotFoundError(f"The file {self.file} does not exist.")
        else:
            self.virus_host = pd.read_csv(self.file, sep='\t')
            self.predictions = self.virus_host[self.virus_host['Predictions'] == 1]
            self.error_check()

    def error_check(self):
        """ Check if the input file is in the correct format. """
        # Check if the file has the three columns needed for visualization
        required_columns = ['pairs', 'InfProbabilities', 'Predictions']
        for column in required_columns:
            if column not in self.virus_host.columns:
                raise ValueError(f"The input file is missing the required column: {column}")
        # Check if the pairs column is in the correct format
        if not all(self.virus_host['pairs'].str.contains(':')):
            raise ValueError("The pairs column is not in the correct format. It should be in the format 'virus:host'")

    def get_unique_virus_host(self):
        """Find the unique virus names and host names in the dataset. The first column contains a virus name and host name, which are separated by a colon."""
        self.unique_viruses = self.virus_host['pairs'].str.split(':').str[0].unique()
        self.unique_hosts = self.virus_host['pairs'].str.split(':').str[1].unique()

    def save_matrix(self):
        """ Save the matrix to a csv file. """
        np.savetxt('virus_host.csv', self.virus_host_array, delimiter=',')
        print('virus_host.csv saved')

    def create_edge_list(self):
        """ Create an edge list from the predictions column of the dataset."""""
        edge_list = []
        for _, row in self.predictions.iterrows():
            virus, host = row['pairs'].split(':')
            edge_list.append((str(virus), str(host)))
        print("Edge list created with {} edges.".format(len(edge_list)))
        # save edge list to a csv file
        edge_df = pd.DataFrame(edge_list, columns=['Virus', 'Host'])
        edge_df.to_csv('edge_list.csv', index=False)
        return edge_list
    
    def graph_from_edge_list(self):
        """ Create a graph from the edge list and draw it using networkx. """
        edge_list = self.create_edge_list()
        G = nx.Graph()
        G.add_edges_from(edge_list)
        plt.figure(figsize=(10, 8))
        pos = nx.spring_layout(G, seed=42)
        nx.draw(G, pos, with_labels=False, node_size=100, node_color= "lightblue" , font_size=10, font_weight='bold', edge_color='gray')
        return G
    
    def graph_igraph(self):
        edge_list = self.create_edge_list()
        # graph the edges
        g = ig.Graph.TupleList(edge_list, directed=False)
        g.vs["label"] = g.vs.indices  # Assign labels to vertices
        layout = g.layout("fr")  # Fruchterman-Reingold layout
        ig.plot(g, layout=layout, vertex_size=300, vertex_color="lightblue", edge_color="gray", vertex_label=g.vs["label"], bbox=(800, 800), margin=50)
        plt.title("Graph from Edge List using igraph")
        plt.show()

    def graph_rustwork(self):
        edge_list = self.create_edge_list()
        # add integer 1 to each edge to indicate presence (for visualization purposes)
        edge_list = [(u, v, int(1)) for u, v in edge_list]  # Add a weight of 1 to each edge
        g = rx.PyGraph()
        g.add_edges_from(edge_list)
        layout = rx.spring_layout(g, seed=42)
        plt.figure(figsize=(10, 8))
        mpl_draw(g)

    
    def centrality_analysis(self, graph):
        print(nx.degree_centrality(graph))
        #find average degree centrality
        degree_centrality = nx.degree_centrality(graph)
        avg_degree_centrality = sum(degree_centrality.values()) / len(degree_centrality)
        print(f"Average Degree Centrality: {avg_degree_centrality:.4f}")


    def initialize_matrix(self, matrix_type:str):
        """"Create a matrix full of zeros and make a list of rows and column labels.
        
        Args:
        matrix_type (str): Type of matrix to create, which determines if the values are boolean or float.
        
        """""
        self.rows = self.unique_viruses
        self.columns = self.unique_hosts
        self.virus_host_array = np.zeros((len(self.rows), len(self.columns)), dtype= bool if matrix_type == "prediction" else float)
    
    def fill_matrix(self, matrix_type:str):
        """ Fill the prediction matrix with 1s and 0s or fill the probability matrix with the InfProbabilities from the dataset."""
        if matrix_type == "prediction":
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
                # Fill the matrix with InfProbabilities from the dataset where virus-host pairs are present
                self.virus_host_array[virus_index][host_index] = row['InfProbabilities']
                #print(row['InfProbabilities'])


    def sort_rows_cols(self, axis:int):
        """ Find the sum of each row and move rows with the highest sum to the top of the matrix. """
        counts = np.sum(self.virus_host_array, axis=axis)
        # Sort inidices based on the counts
        sorted_indices = np.argsort(counts)[::-1]
        # Sort the matrix and row names based on the sorted indices
        if axis == 1:
            self.rows = self.rows[sorted_indices]
            self.virus_host_array = self.virus_host_array[sorted_indices]
        elif axis == 0:
            self.columns = self.columns[sorted_indices]
            self.virus_host_array = self.virus_host_array[:, sorted_indices]


    def sort_matrix(self):
        self.sort_rows_cols(1)
        self.sort_rows_cols(0)

    def expand_matrix(self, matrix_type:str):
        """ Expand the matrix to make it square by adding rows for the hosts and columns for the viruses.
        
        Args:
        matrix_type (str): Type of matrix to create, which determines if the values are boolean or float.
        
        """
        # add number of hosts to the rows
        self.virus_host_array_square = np.concatenate((self.virus_host_array, np.zeros((len(self.unique_hosts), len(self.columns)), dtype=bool if matrix_type == "prediction" else float)), axis=0)
        self.rows_square = np.concatenate((self.rows, self.unique_hosts))
        # add number of viruses to the columns
        self.virus_host_array_square = np.concatenate((self.virus_host_array_square, np.zeros((len(self.rows_square), len(self.unique_viruses)), bool if matrix_type == "prediction" else float)), axis=1)
        self.columns_square = np.concatenate((self.columns, self.unique_viruses))
        # Fill with 1s
        self.fill_bottom_right(matrix_type)

    def fill_bottom_right(self, matrix_type:str):
        """ Fill the virus-host pairs with same values as host-virus pairs. Use the predictions or the InfProbabilities from the dataset."""

        if matrix_type == "prediction":
            arrayname = self.predictions
        else:
            arrayname = self.virus_host

        for _, row in arrayname.iterrows():
            virus, host = row['pairs'].split(':')
            virus_index = np.where(self.columns_square == virus)[0][0]
            host_index = np.where(self.rows_square == host)[0][0]
            self.virus_host_array_square[host_index][virus_index] = 1 if matrix_type == "prediction" else row['InfProbabilities']
        return self.virus_host_array_square
    
    def make_rectangular_matrix(self, matrix_type:str):
        """ Call all functions to make a rectangular matrix"""
        matrix_type1 = matrix_type.lower()
        self.get_unique_virus_host()
        self.initialize_matrix(matrix_type1)
        self.fill_matrix(matrix_type1)
        self.sort_matrix()
        return self.virus_host_array

    def make_square_matrix(self, matrix_type:str):
        """ Call all functions to make a square matrix"""
        self.make_rectangular_matrix(matrix_type.lower())
        self.expand_matrix(matrix_type.lower())
        return self.virus_host_array_square
    
    # 576 partner task
    def square_to_rectangular_matrix(self, matric_type:str):
        """Convert the Square matrix back into the original rectangular matrix"""
        self.get_unique_virus_host()
        unique_virus_count = len(set(self.unique_viruses))
        unique_host_count = len(set(self.unique_hosts))
        self.virus_host_array = self.virus_host_array[:unique_virus_count, :unique_host_count]
        return self.virus_host_array

    def plot_heatmap(self, matrix_type:str):
        # if matrix type prediction, use purples, otherwise use warm colors
        # Make heatmap color red if 1, grey if 0.5, and blue if 0
        # Make a user-defined colormap.
        cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",["r", "white", "b"])
        sns.heatmap(self.virus_host_array, cmap= mcol.LinearSegmentedColormap.from_list("MyCmapName",["white", "cadetblue"]) if matrix_type == 'prediction' else cm1)
        plt.gcf().set_size_inches(7, 14)
        plt.xlabel("Hosts")
        plt.ylabel("Viruses")
        matrix_title = self.file.replace('Sample_Input/', '')
        plt.title(matrix_title.split('_')[0] + ' ' + matrix_type)
        # save the figure in the heatmaps folder
        # get matrix name before the first underscore
        plt.savefig('Heatmaps/Heatmap_' + matrix_title.split('_')[0] + '_' +  matrix_type + '.png')
        plt.show()

