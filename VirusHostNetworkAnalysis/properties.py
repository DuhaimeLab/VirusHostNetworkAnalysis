import networkx as nx
import matplotlib.pyplot as plt
from multiprocessing import Pool
import numpy as np
import matplotlib.colors as mcol
import seaborn as sns
from tqdm import tqdm
from networkx.algorithms import community
from typing import List, Tuple
from VirusHostNetworkAnalysis.prediction_matrix import PredictionMatrix


# change to properties.py
# in the commments include what the arguments are
# change name of graph to BipartiteGraph
# make the degree distribution a line
# is there a way to return the plot? if you return it then you can add stuff on top of it


class BipartiteGraph:
    """ Take in a matrix and create a graph from it. The graph is initialized with a given number of rows and columns.

    Args:
    MatrixClass: Can be a PredictionMatrix, ER, or ConfigurationModel class.
    
    """
    def __init__(self, MatrixClass):
        self.MatrixClass = MatrixClass
        self.input_matrix = np.array(MatrixClass.virus_host_array)
        self.rows = MatrixClass.rows
        self.columns = MatrixClass.columns
        if isinstance(MatrixClass, PredictionMatrix):
            self.predictions = MatrixClass.predictions
            self.probability = MatrixClass.probability
            self.title = MatrixClass.title
        else: 
            # Most of the time the input matrix is a boolean matrix. For example, ER and Configuration models will always be boolean.
            self.probability = False
            self.title = "Null_Model"
    
    def sort_rows_cols(self, axis:int):
        """ Find the sum of each row and move rows with the highest sum to the top of the matrix. 
        
        Args:
        axis (int): Axis to be used when sorting the matrix. 0 for rows, 1 for columns.
        
        """
        # Calculate the sum of each row or column
        counts = np.sum(self.input_matrix, axis=axis)
        # Sort inidices based on the counts. Sorts the sums in descending order, so any ties will be broken by descending order of the index.
        sorted_indices = np.argsort(counts)[::-1]
        # Sort the matrix and row names based on the sorted indices
        if axis == 1:
            self.rows = self.rows[sorted_indices]
            self.input_matrix = self.input_matrix[sorted_indices]
        elif axis == 0:
            self.columns = self.columns[sorted_indices]
            self.input_matrix = self.input_matrix[:, sorted_indices]

    def sort_matrix(self):
        """ Sort the matrix by rows and columns. """
        # Sort the rows
        self.sort_rows_cols(1)
        # Sort the columns
        self.sort_rows_cols(0)


    def fill_bottom_right(self):
        """ Fill the virus-host pairs with same values as host-virus pairs. 
        Use the predictions or the InfProbabilities from the dataset.
        """
        if self.probability is False:
            arrayname = self.predictions
        else:
            arrayname = self.MatrixClass.virus_host

        for _, row in arrayname.iterrows():
            virus, host = row['pairs'].split(':')
            virus_index = np.where(self.columns_square == virus)[0][0]
            host_index = np.where(self.rows_square == host)[0][0]
            self.virus_host_array_square[host_index][virus_index] = 1 if self.probability is False else row['InfProbabilities'or'Scores']
        return self.virus_host_array_square
    
    def make_square_matrix(self):
        """ Expand the matrix to make it square by adding rows for the hosts and columns for the viruses.
        """
        # rectangular matrix should be sorted first
        self.sort_matrix()
        # add number of hosts to the rows
        self.virus_host_array_square = np.concatenate((self.input_matrix, np.zeros((len(self.columns), len(self.columns)), dtype=bool if self.probability is False else float)), axis=0)
        self.rows_square = np.concatenate((self.rows, self.columns))
        # add number of viruses to the columns
        self.virus_host_array_square = np.concatenate((self.virus_host_array_square, np.zeros((len(self.rows_square), len(self.rows)), bool if self.probability is False else float)), axis=1)
        self.columns_square = np.concatenate((self.columns, self.rows))
        # Fill with 1s
        self.fill_bottom_right()
        return self.virus_host_array_square

    
    def initialize_graph(self):
        """ Initialize the graph with nodes and edges. """
        self.G = nx.Graph()
        self.G.add_nodes_from(self.rows)
        self.G.add_nodes_from(self.columns)
        # self.G.add_edges_from(np.argwhere(self.input_matrix == 1))
        # Add edges between nodes based on the input matrix
        for i in range(len(self.input_matrix)):
            for j in range(len(self.input_matrix[1])):
                if self.input_matrix[i][j] == 1:
                    self.G.add_edge(self.rows[i], self.columns[j])

    def draw_graph(self, include_label:bool):
        """ Draw the graph using NetworkX. 
        Args:
        include_label (bool): Whether to include labels on the nodes or not.
        """

        plt.figure(figsize=(20, 15))
        self.initialize_graph()
        # Draw the graph. Blue for nodes in rows, red for nodes in columns
        node_color = ['blue' if node in self.rows else 'red' for node in self.G.nodes()]

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
            self.eigenvector_virus = {k: v for k, v in self.eigenvector.items() if k in self.rows}
            self.eigenvector_host = {k: v for k, v in self.eigenvector.items() if k in self.columns}
            print("eigen done")

        elif algorithm == "betweenness":
            self.betweenness = nx.betweenness_centrality(self.G)
            self.betweenness_virus = {k: v for k, v in self.betweenness.items() if k in self.rows}
            self.betweenness_host = {k: v for k, v in self.betweenness.items() if k in self.columns}
            print("betweenness done")

        elif algorithm == "closeness":
            self.closeness = nx.closeness_centrality(self.G)
            self.closeness_virus = {k: v for k, v in self.closeness.items() if k in self.rows}
            self.closeness_host = {k: v for k, v in self.closeness.items() if k in self.columns}
            print("closeness done")
        
        else:
            raise ValueError("Algorithm not supported. Choose from 'eigenvector', 'betweenness', or 'closeness'.")
    

    def plot_degree_distribution(self):
        """ Plot the degree distribution of the graph. """
        degree_seq = self.calculate_degree()

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
       
    # change name to plot_eigenvector_centrality 
    def plot_eigenvector_centrality(self):
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

    def plot_betweenness_centrality(self):
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

    def plot_closeness_centrality(self):
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

        # Sort the matrix so that patterns are visible
        self.sort_matrix()

        # Make heatmap color red if 1, grey if 0.5, and blue if 0 using user-defined color map
        cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",[(ranges[0], color_map[0]), (ranges[1], color_map[1]), (ranges[2], color_map[2]), (ranges[3], color_map[2]), (ranges[4], color_map[3]), (ranges[5], color_map[4])])
        sns.heatmap(self.input_matrix, cmap= mcol.LinearSegmentedColormap.from_list("MyCmapName",["white", prediction_color]) if self.probability is False else cm1)
        plt.gcf().set_size_inches(7, 14)
        plt.xlabel("Hosts")
        plt.ylabel("Viruses")
        plt.title(self.title)
        # save the figure in the heatmaps folder
        # get matrix name before the first underscore
        plt.savefig('Heatmaps/Heatmap_' + self.title + '_' + 'predictions' if self.probability is False else 'probabiltiies' + '.png')
        plt.show()


### NESTEDNESS FIX 
    def nestedness_rows(self, pair):
        """Calculate nestedness for rows using the NODF algorithm for the array."""
        # generate list of rows to compare
        pair1 = self.input_matrix[pair[0],]
        pair2 = self.input_matrix[pair[1],]
        N_row = self.compare(pair1, pair2)
        #print(pair[0], pair[1], N_row)
        #print(N_row)
        return N_row

    def nestedness_cols(self, pair):
        """Calculate nestedness for cols using the NODF algorithm for the array."""
        pair1 = self.input_matrix[:, pair[0]]
        pair2 = self.input_matrix[:, pair[1]]
        N_col = self.compare(pair1, pair2)  # pyright: ignore
        return N_col

    def pairs(self, axis: int = 0) -> List[Tuple[int, int]]:
        """Determine all possible i-j pairs.

        Args:
            axis (int): Axis to be used when determining all pairs.
        """
        lst: List[Tuple[int, int]] = []
        for i in range(0, self.input_matrix.shape[axis]):
            for j in range(i + 1, self.input_matrix.shape[axis]):
                lst.append((i, j))
        return lst


    def compare(self, x: list[int], y: list[int]) -> float:
        """Compare two lists containing 0 and 1.

        Args:
            x (list[int]): first list
            y (list[int]): second list
        """
        if sum(x) <= sum(y):
            val = 0
        elif sum(y) == 0:
            val = 0
        else:
            counter = 0
            total = 0
            for i, j in zip(x, y):
                if i == 1 and j == 1:
                    counter += 1
                    total += 1
                elif i == 0 and j == 1:
                    total += 1
            val = counter / total
        return val * 100

    def run_parallel(self, num_procs: int = 6):
        """Run multiple process of the compute_feature method.

        This significantly improves run time by using multiple CPU cores.

        Args:
            num_procs (int): Number of core to be used.
        """

        # Sort the matrix
        self.sort_matrix()

        with Pool(num_procs) as pool:
            #nrow = pool.map(self.nestedness_rows, self.pairs(axis=0))
            with tqdm(total=len(self.pairs(axis=0)) + len(self.pairs(axis=1)), desc="Calculating nestedness", colour="green") as pbar:
                nrow = []
                for result in pool.imap_unordered(self.nestedness_rows, self.pairs(axis=0)):
                    nrow.append(result)
                    pbar.update(1)

                ncol = []
                for result in pool.imap_unordered(self.nestedness_cols, self.pairs(axis=1)):
                    ncol.append(result)
                    pbar.update(1)

        nrow = sum(nrow) / (len(self.input_matrix) * (len(self.input_matrix) - 1) / 2)
        ncol = sum(ncol) / (len(self.input_matrix[0]) * (len(self.input_matrix[0]) - 1) / 2)
        print(nrow, ncol)
        nodf = (nrow + ncol) / 2
        return nodf
    


    def calculate_degree(self):
        """ Calculate the degree of each node in the graph. """
        host_degrees = []
        for col in range(len(self.input_matrix[0])):
            total = 0
            for row in range(len(self.input_matrix)):
                if self.input_matrix[row][col] == 1:
                    total += 1
            host_degrees.append(total)
        
        virus_degrees = []
        for row in range(len(self.input_matrix)):
            total = 0
            for col in range(len(self.input_matrix[0])):
                if self.input_matrix[row][col] == 1:
                    total += 1
            virus_degrees.append(total)

        return [virus_degrees, host_degrees]
    

    # Calculate modularity
    def calculate_modularity(self):
        """ Calculate the modularity of the graph. """
        G = nx.from_numpy_array(self.virus_host_array_square)
        partition = community.greedy_modularity_communities(G)
        modularity = community.modularity(G, partition)
        print("Modularity: ", modularity)
        return modularity