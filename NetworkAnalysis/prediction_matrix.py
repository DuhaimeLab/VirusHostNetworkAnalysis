import numpy as np
import pandas as pd
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt

## Square Matrix Class
class PredictionMatrix:
    def __init__(self, vhip_predictions_file:str):
        self.file = vhip_predictions_file
        self.virus_host = pd.read_csv(self.file, sep='\t')
        self.predictions = self.virus_host[self.virus_host['Predictions'] == 1]

    # Find the unique virus names and host names in the dataset. The first column contains a virus name and host name, which are separated by a colon.
    def get_unique_virus_host(self):
        self.unique_viruses = self.virus_host['pairs'].str.split(':').str[0].unique()
        self.unique_hosts = self.virus_host['pairs'].str.split(':').str[1].unique()
        print('Unique viruses and hosts found')

    # Save the matrix to a csv file
    def save_matrix(self):
        np.savetxt('virus_host.csv', self.virus_host_array, delimiter=',')
        print('virus_host.csv saved')

    # Create a matrix full of zeros and make a list of rows and column labels
    def initialize_matrix(self):
        self.rows = self.unique_viruses
        self.columns = self.unique_hosts
        self.virus_host_array = np.zeros((len(self.rows), len(self.columns)), dtype=bool)
        print('Matrix initialized') 
    
    def fill_matrix(self):
        # Loop through the rows in subset of data where predictions is 1
        for _, row in self.predictions.iterrows():
            virus, host = row['pairs'].split(':')
            # Find the index of the virus (row) and host (col) in the matrix
            virus_index = np.where(self.rows == virus)[0][0]
            host_index = np.where(self.columns == host)[0][0]
            # Fill the matrix with 1
            self.virus_host_array[virus_index][host_index] = 1
        print('Matrix filled')

    def sort_rows(self):
        # Count the 1s in each row
        row_counts = np.sum(self.virus_host_array, axis=1)
        # Sort inidices based on the counts
        sorted_indices = np.argsort(row_counts)[::-1]
        # Sort the matrix and row names based on the sorted indices
        self.virus_host_array = self.virus_host_array[sorted_indices]
        self.rows = self.rows[sorted_indices]
        print('rows sorted')

    def sort_columns(self):
        # Count the 1s in each column
        col_counts = np.sum(self.virus_host_array, axis=0)
        # Sort inidices based on the counts
        sorted_indices = np.argsort(col_counts)[::-1]
        # Sort the matrix and column names based on the sorted indices
        self.virus_host_array = self.virus_host_array[:, sorted_indices]
        self.columns = self.columns[sorted_indices]
        print('columns sorted')

    def sort_matrix(self):
        self.sort_rows()
        self.sort_columns()

    # Make the matrix a square by adding rows for the hosts and columns for the viruses
    def expand_matrix(self):
        # add number of hosts to the rows
        self.virus_host_array_square = np.concatenate((self.virus_host_array, np.zeros((len(self.unique_hosts), len(self.columns)), dtype=bool)), axis=0)
        self.rows_square = np.concatenate((self.rows, self.unique_hosts))
        # add number of viruses to the columns
        self.virus_host_array_square = np.concatenate((self.virus_host_array_square, np.zeros((len(self.rows_square), len(self.unique_viruses)), dtype=bool)), axis=1)
        self.columns_square = np.concatenate((self.columns, self.unique_viruses))
        # Fill with 1s
        self.fill_bottom_right()

    def fill_bottom_right(self):
        # fill the virus-host pairs with same values as host-virus pairs
        for index, row in self.predictions.iterrows():
            virus, host = row['pairs'].split(':')
            virus_index = np.where(self.columns_square == virus)[0][0]
            host_index = np.where(self.rows_square == host)[0][0]
            self.virus_host_array_square[host_index][virus_index] = 1
        return self.virus_host_array_square
    
    # Make rectangular matrix
    def make_rectangular_matrix(self):
        self.get_unique_virus_host()
        self.initialize_matrix()
        self.fill_matrix()
        self.sort_matrix()
        return self.virus_host_array

    # Make square matrix
    def make_square_matrix(self):
        self.make_rectangular_matrix()
        self.expand_matrix()
        return self.virus_host_array_square

    # Plot the heatmap for the rectangular matrix
    def plot_heatmap(self):
        sns.heatmap(self.virus_host_array, cmap="Purples")
        # change height of heatmap
        plt.gcf().set_size_inches(7, 14)
        plt.xlabel("Hosts")
        plt.ylabel("Viruses")
        plt.savefig("heatmap.png")
        plt.show()
