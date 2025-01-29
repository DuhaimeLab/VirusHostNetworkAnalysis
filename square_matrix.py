import numpy as np
import pandas as pd
# # read in the tsv file
# virus_host = pd.read_csv('Aug4_predictions.tsv', sep='\t')

# # The first column contains a virus name and host name, which are separated by a colon. 
# # Find the unique virus names and host names in the dataset.
# unique_viruses = virus_host['pairs'].str.split(':').str[0].unique()
# unique_hosts = virus_host['pairs'].str.split(':').str[1].unique()

# print(len(unique_viruses))
# print(len(unique_hosts))

# virus_host_array = np.zeros((len(unique_viruses)+len(unique_hosts), len(unique_hosts)+len(unique_viruses)), dtype=bool)
# rows = np.concatenate((unique_viruses, unique_hosts))
# columns = np.concatenate((unique_hosts, unique_viruses))
# print("check")

# # Subset of data where predictions is 1
# subset_virus_host = virus_host[virus_host['Predictions'] == 1]
# # Loop through the rows in subset
# for index, row in subset_virus_host.iterrows():
#     virus, host = row['pairs'].split(':')
#     # Find the index of the virus (row) and host (col) in the matrix
#     virus_index = np.where(rows == virus)[0][0]
#     host_index = np.where(columns == host)[0][0]
#     # Fill the matrix with 1
#     virus_host_array[virus_index][host_index] = 1

#     # Find the index of the virus (col) and host (row) in the matrix
#     virus_index = np.where(columns == virus)[0][0]
#     host_index = np.where(rows == host)[0][0]
#     # Fill with 1
#     virus_host_array[host_index][virus_index] = 1

# print("done with loop")
# print(virus_host_array)

# # Save to csv
# np.savetxt('virus_host.csv', virus_host_array, delimiter=',')
# print("saved")

# # import networkx as nx
# # import matplotlib.pyplot as plt
# # G = nx.Graph(virus_host_array)
# # plt.figure(figsize=(12, 12))
# # pos = nx.random_layout(G)
# # # make light blue if virus, blue if host
# # nx.draw(G, pos, with_labels=False, node_size=50, node_color="lightblue", edge_color='gray')
# # #remove any nodes that are not connected to any other nodes
# # G.remove_nodes_from(list(nx.isolates(G)))
# # plt.show()

# # Create empty graph with nodes for virus and host
# G = np.Graph()
# for virus in unique_viruses:
#     G.add_node(virus, type='virus')
# for host in unique_hosts:
#     G.add_node(host, type='host')

# # ER graph
# for row in rows:
#     for col in columns:
#         if row != col:
#             # set p to a random value between 0 and 1
#             p = np.random.rand()
#             # if p is less than 0.1, add an edge between the nodes
#             if p < 0.1:
#                 G.add_edge(row, col)
                

class square_matrix:
    def __init__(self, vhip_predictions_file:str):
        self.file = vhip_predictions_file
        self.virus_host = pd.read_csv(self.file, sep='\t')
        self.predictions = self.virus_host[self.virus_host['Predictions'] == 1]

    # Find the unique virus names and host names in the dataset. 
    # The first column contains a virus name and host name, which are separated by a colon.
    def get_unique_virus_host(self):
        self.unique_viruses = self.virus_host['pairs'].str.split(':').str[0].unique()
        self.unique_hosts = self.virus_host['pairs'].str.split(':').str[1].unique()
    
    # Create a matrix full of zeros and make a list of rows and column labels
    def initialize_matrix(self):
        self.rows = np.concatenate((self.unique_viruses, self.unique_hosts))
        self.columns = np.concatenate((self.unique_hosts, self.unique_viruses))
        self.virus_host_array = np.zeros((len(self.rows), len(self.columns)), dtype=bool)

    # Fill the matrix with the predictions values from the dataset.
    def fill_matrix(self):
        # Loop through the rows in subset of data where predictions is 1
        for index, row in self.predictions.iterrows():
            virus, host = row['pairs'].split(':')
            # Find the index of the virus (row) and host (col) in the matrix
            virus_index = np.where(self.rows == virus)[0][0]
            host_index = np.where(self.columns == host)[0][0]
            # Fill the matrix with 1
            self.virus_host_array[virus_index][host_index] = 1
            # Find the index of the virus (col) and host (row) in the matrix
            virus_index = np.where(self.columns == virus)[0][0]
            host_index = np.where(self.rows == host)[0][0]
            # Fill with 1
            self.virus_host_array[host_index][virus_index] = 1

    # Save the matrix to a csv file
    def save_matrix(self):
        np.savetxt('virus_host.csv', self.virus_host_array, delimiter=',')

    # Call all the funcitons above to make the matrix
    def make_square_matrix(self):
        self.get_unique_virus_host()
        self.initialize_matrix()
        self.fill_matrix()
        self.save_matrix()

