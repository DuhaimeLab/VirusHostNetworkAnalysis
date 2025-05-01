from calendar import c
import networkx as nx
import matplotlib.pyplot as plt
from multiprocessing import Pool
from functools import partial
import numpy as np
import matplotlib.colors as mcol
import seaborn as sns
from multiprocessing import Pool
import itertools

class Graph:
    """ Take in a matrix and create a graph from it. The graph is initialized with a given number of rows and columns.

    Args:
    input_matrix (np.ndarray): The matrix to create the graph from.
    
    """
    def __init__(self, input_matrix, virus_labels, host_labels):
        self.input_matrix = input_matrix
        self.x_labels = virus_labels
        self.y_labels = host_labels
        
        #print(type(self.input_matrix))
    
    def initialize_graph(self):
        """ Initialize the graph with nodes and edges. """
        self.G = nx.Graph()
        self.G.add_nodes_from(self.x_labels)
        self.G.add_nodes_from(self.y_labels)
        # self.G.add_edges_from(np.argwhere(self.input_matrix == 1))
        # Add edges between nodes based on the input matrix
        for i in range(len(self.input_matrix)):
            for j in range(len(self.input_matrix[1])):
                if self.input_matrix[i][j] == 1:
                    self.G.add_edge(self.x_labels[i], self.y_labels[j])


    def draw_graph(self, include_label:bool):
        """ Draw the graph using NetworkX. """
        plt.figure(figsize=(20, 15))
        self.initialize_graph()
        # Draw the graph. Blue for nodes in rows, red for nodes in columns
        node_color = ['blue' if node in self.x_labels else 'red' for node in self.G.nodes()]

        # Add one to that nodes with 0 degrees are visible
        #node_size = [(self.G.degree(node)+1) * 50 for node in self.G.nodes()]

        # Set node size proportional to the degree of the node
        pos = nx.random_layout(self.G, seed=42)
        nx.draw(self.G, pos, with_labels= include_label, node_color=node_color,
                node_size = 100) 
        
    
    # Calculate the centrality of the graph
    def calculate_centrality(self, max_iter, algorithm = "eigenvector"):
        """ Calculate the centrality of the graph. """
        if algorithm == "eigenvector":
            self.eigenvector = nx.eigenvector_centrality(self.G, max_iter=max_iter)
            # calculate eigenvector centrality for only the virus and only host
            self.eigenvector_virus = {k: v for k, v in self.eigenvector.items() if k in self.x_labels}
            self.eigenvector_host = {k: v for k, v in self.eigenvector.items() if k in self.y_labels}
            print("eigen done")

        elif algorithm == "betweenness":
            self.betweenness = nx.betweenness_centrality(self.G)
            self.betweenness_virus = {k: v for k, v in self.betweenness.items() if k in self.x_labels}
            self.betweenness_host = {k: v for k, v in self.betweenness.items() if k in self.y_labels}
            print("betweenness done")

        elif algorithm == "closeness":
            self.closeness = nx.closeness_centrality(self.G)
            self.closeness_virus = {k: v for k, v in self.closeness.items() if k in self.x_labels}
            self.closeness_host = {k: v for k, v in self.closeness.items() if k in self.y_labels}
            print("closeness done")
        
        else:
            raise ValueError("Algorithm not supported. Choose from 'eigenvector', 'betweenness', or 'closeness'.")
    

    def degree_distribution(self, degree_seq):
        """ Plot the degree distribution of the graph. """
        plt.figure(figsize=(10, 6))
        virus_degree = degree_seq[0]
        host_degree = degree_seq[1]

        # Plot histograms for virus and host degree distributions side by side)
        plt.hist(virus_degree, color='g', density=True)
        plt.title('Virus Degree Distribution')
        plt.xlabel('Degree')
        plt.ylabel('Frequency')
        plt.grid()
        plt.show()
    
        plt.hist(host_degree, color='b', density=True)
        plt.title('Host Degree Distribution')
        plt.xlabel('Degree')
        plt.ylabel('Frequency')
        plt.grid()
        plt.show()
       
    def plot_eigenvectors(self):
        """ Plot the eigenvector centrality of the graph. """
      
        # Plot for the virus eigenvector centrality
        plt.figure(figsize=(10, 6))
        plt.hist(list(self.eigenvector_virus.values()), color='g', density=True)
        plt.title('Eigenvector Centrality for Viruses')
        plt.xlabel('Eigenvector Centrality')
        plt.ylabel('Frequency')
        plt.grid()
        plt.show()

        # Plot for the host eigenvector centrality
        plt.figure(figsize=(10, 6))
        plt.hist(list(self.eigenvector_host.values()), color='b', density=True)
        plt.title('Eigenvector Centrality for Hosts')
        plt.xlabel('Eigenvector Centrality')
        plt.ylabel('Frequency')
        plt.grid()
        plt.show()

    def plot_betweenness(self):
        """ Plot the betweenness centrality of the graph. """

        # Plot for the virus betweenness centrality
        plt.figure(figsize=(10, 6))
        plt.hist(list(self.betweenness_virus.values()), color='g', density=True)
        plt.title('Betweenness Centrality for Viruses')
        plt.xlabel('Betweenness Centrality')
        plt.ylabel('Frequency')
        plt.grid()
        plt.show()

        # Plot for the host betweenness centrality
        plt.figure(figsize=(10, 6))
        plt.hist(list(self.betweenness_host.values()), color='b', density=True)
        plt.title('Betweenness Centrality for Hosts')
        plt.xlabel('Betweenness Centrality')
        plt.ylabel('Frequency')
        plt.grid()
        plt.show()

    def plot_closeness(self):
        """ Plot the closeness centrality of the graph. """

        # Plot for the virus closeness centrality
        plt.figure(figsize=(10, 6))
        plt.hist(list(self.closeness_virus.values()), color='g', density=True)
        plt.title('Closeness Centrality for Viruses')
        plt.xlabel('Closeness Centrality')
        plt.ylabel('Frequency')
        plt.grid()
        plt.show()

        # Plot for the host closeness centrality
        plt.figure(figsize=(10, 6))
        plt.hist(list(self.closeness_host.values()), color='b', density=True)
        plt.title('Closeness Centrality for Hosts')
        plt.xlabel('Closeness Centrality')
        plt.ylabel('Frequency')
        plt.grid()
        plt.show()

    def closeness_boxplot(self, iterations):
        # iterations is a list of dictionaries with the closeness centrality for each iteration
        closeness_virus = []
        closeness_host = []
        for i in range(iterations):
            closeness_virus.append(list(self.closeness_virus[i].values()))
            closeness_host.append(list(self.closeness_host[i].values()))
        plt.figure(figsize=(10, 6))
        
    def plot_heatmap(self, prediction_color = "indigo", color_map=["red", "lightpink", "white", "#a2cffe", "blue"], ranges=[0, 0.2, 0.45, 0.55, 0.8, 1]):
        """ Plot the heatmap of the matrix.
        Let the user choose the colors and ranges of the heatmap.
        """

        # Matrix type is 'prediction' if matrix only has 1s and 0s
        matrix_type = 'prediction' if np.all(np.isin(self.input_matrix, [0, 1])) else 'probability'

        # Make heatmap color red if 1, grey if 0.5, and blue if 0 using user-defined color map
        cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",[(ranges[0], color_map[0]), (ranges[1], color_map[1]), (ranges[2], color_map[2]), (ranges[3], color_map[2]), (ranges[4], color_map[3]), (ranges[5], color_map[4])])
        sns.heatmap(self.input_matrix, cmap= mcol.LinearSegmentedColormap.from_list("MyCmapName",["white", prediction_color]) if matrix_type == 'prediction' else cm1)
        plt.gcf().set_size_inches(7, 14)
        plt.xlabel("Hosts")
        plt.ylabel("Viruses")
        matrix_title = "Config Model"
        plt.title(matrix_title.split('_')[0] + ' ' + matrix_type)
        # save the figure in the heatmaps folder
        # get matrix name before the first underscore
        plt.savefig('Heatmaps/Heatmap_' + matrix_title.split('_')[0] + '_' +  matrix_type + '.png')
        plt.show()