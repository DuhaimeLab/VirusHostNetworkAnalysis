import matplotlib
import networkx as nx
import matplotlib.pyplot as plt
from multiprocessing import Pool
import numpy as np
import matplotlib.colors as mcol
import seaborn as sns
from typing import List, Tuple
from VirusHostNetworkAnalysis.prediction_matrix import PredictionMatrix
import infomap

# in the commments include what the arguments are
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
            #self.predictions = MatrixClass.predictions
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
            arrayname = self.MatrixClass.virus_host[self.MatrixClass.virus_host['Predictions'] == 1]
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

        # Set node size proportional to the degree of the node
        pos = nx.random_layout(self.G, seed=42)
        nx.draw(self.G, pos, with_labels= include_label, node_color=node_color,
                node_size = 100) 
        
    def unipartite_matrix(self):
        """ Create a unipartite matrix from the bipartite matrix. 
        The matrix needs to be sorted for the algorithm to work."""
        # Viruses
        self.unipartite_viruses = np.zeros((len(self.rows), len(self.rows)), dtype=int)
        self.unipartite_hosts = np.zeros((len(self.columns), len(self.columns)), dtype=int)
        for col in self.input_matrix.T:
            # to save time, if the sum of the row is 0 or 1, skip it
            if sum(col) > 1:
                for i in range(len(col)):
                # if spot is 0, do not need to check the rest of the column
                    if col[i] == 1:
                        for j in range(i+1, len(col)):
                            if col[j] == 1:
                                self.unipartite_viruses[i][j] += 1
                                self.unipartite_viruses[j][i] += 1
            else:
                #end the loop if the sum of the column is 0 or 1
                # assuming that the matrix is sorted, this means that the rest of the column will be 0
                break
        # Hosts
        for row in self.input_matrix:
            # to save time, if the sum of the row is 0 or 1, skip it
            if sum(row) > 1:
                for i in range(len(row)):
                    # if spot is 0, do not need to check the rest of the row
                    if row[i] == 1:
                        for j in range(i+1, len(row)):
                            if row[j] == 1:
                                # only fill to the upper triangle of the matrix
                                self.unipartite_hosts[i][j] += 1
                                self.unipartite_hosts[j][i] += 1
            else:
                # end the loop if the sum of the row is 0 or 1
                # assuming that the matrix is sorted, this means that the rest of the row will be 0
                break

    def unipartite_graph(self):
        # # Make a graph from the unipartite matrix for hosts only
        self.G_host = nx.Graph()
        self.G_host.add_nodes_from(self.columns)
        for i in range(len(self.unipartite_hosts)):
            for j in range(i+1, len(self.unipartite_hosts[1])):
                if self.unipartite_hosts[i][j] > 0:
                    # scale weight so that the maximum weight is 1 aka divide by the maximum weight
                    self.G_host.add_edge(self.columns[i], self.columns[j], weight=self.unipartite_hosts[i][j]/np.max(self.unipartite_hosts))
        plt.figure(figsize=(20, 20))
        pos_host = nx.circular_layout(self.G_host, scale=2)
        # add an alpha to the edges so that they are not too dark
        nx.draw(self.G_host, pos_host, with_labels=False, node_color='blue', node_size=400, width=[4*np.power(self.G_host[u][v]['weight'], 2) for u, v in self.G_host.edges()], alpha =0.6)
                
    # Calculate the centrality of the graph
    def calculate_centrality(self, algorithm = "eigenvector", max_iter = 1000):
        """ Calculate the centrality of the graph. """
        self.initialize_graph()
        if algorithm == "eigenvector":
            self.eigenvector = nx.eigenvector_centrality(self.G, max_iter=max_iter)
            # Split the eigenvector centrality into virus and host
            self.eigenvector_virus = {k: v for k, v in self.eigenvector.items() if k in self.rows}
            self.eigenvector_host = {k: v for k, v in self.eigenvector.items() if k in self.columns}
        elif algorithm == "betweenness":
            self.betweenness = nx.betweenness_centrality(self.G)
            # Split the betweenness centrality into virus and host
            self.betweenness_virus = {k: v for k, v in self.betweenness.items() if k in self.rows}
            self.betweenness_host = {k: v for k, v in self.betweenness.items() if k in self.columns}
        elif algorithm == "closeness":
            self.closeness = nx.closeness_centrality(self.G)
            # Split the closeness centrality into virus and host
            self.closeness_virus = {k: v for k, v in self.closeness.items() if k in self.rows}
            self.closeness_host = {k: v for k, v in self.closeness.items() if k in self.columns}   
        else:
            raise ValueError("Algorithm not supported. Choose from 'eigenvector', 'betweenness', or 'closeness'.")

    def calculate_degree(self):
        """ Calculate the degree of each node in the graph. """
        
        host_degrees = []
        for col in range(len(self.input_matrix[0])):
            total = 0
            for row in range(len(self.input_matrix)):
                # for each column/host, check how many viruses are connected to it
                if self.input_matrix[row][col] == 1:
                    total += 1
            host_degrees.append(total)
        
        virus_degrees = []
        for row in range(len(self.input_matrix)):
            total = 0
            for col in range(len(self.input_matrix[0])):
                # for each row/virus, check how many hosts are connected to it
                if self.input_matrix[row][col] == 1:
                    total += 1
            virus_degrees.append(total)

        return [virus_degrees, host_degrees]        

    def plot_degree_distribution(self):
        """ Plot the degree distribution of the graph. """
        # Calculate the degree of each node in the graph
        degree_seq = self.calculate_degree()
        plt.figure(figsize=(10, 6))
        # Divide into two lists, one for virus and one for host
        virus_degree = degree_seq[0]
        host_degree = degree_seq[1]

        # Plot histograms for virus and host degree distributions
        # Virus distribution is green and host distribution is blue
        sns.kdeplot(virus_degree, color='g', fill=True) 
        plt.title('Virus Degree Distribution')
        plt.xlabel('Degree')
        plt.ylabel('Frequency')
        plt.grid()
        plt.show()

        sns.kdeplot(host_degree, color='b', fill=True)
        plt.title('Host Degree Distribution')
        plt.xlabel('Degree')
        plt.ylabel('Frequency')
        plt.grid()
        plt.show()

    def plot_degree_by_species(self):
        """ Plot the degree of each node in the graph. """
        # Calculate the degree of each node in the graph
        degree_seq = self.calculate_degree()
        plt.figure(figsize=(10, 6))
        # Divide into two lists, one for virus and one for host
        virus_degree = degree_seq[0]
        host_degree = degree_seq[1]

        # dictionary of virus names and their corresponding degree
        virus_dict = {self.rows[i]: virus_degree[i] for i in range(len(virus_degree))}
        # dictionary of host names and their corresponding degree
        host_dict = {self.columns[i]: host_degree[i] for i in range(len(host_degree))}

        # sort the dictionaries by value
        virus_dict = dict(sorted(virus_dict.items(), key=lambda item: item[1], reverse=True))
        host_dict = dict(sorted(host_dict.items(), key=lambda item: item[1], reverse=True))

        # plot virus name on x axis and corresponding degree on y axis
        plt.figure(figsize=(10, 6))
        # dot plot for virus
        plt.scatter(list(virus_dict.keys()), list(virus_dict.values()), color='g')
        # remove x axis labels
        plt.xticks([])
        plt.title('Virus Degree Distribution')
        plt.xlabel('Virus')
        plt.ylabel('Degree')
        # Add Generalist and Specialist labels
        plt.text(-0.5, -2, 'Generalist', fontsize=12, ha='right', va='center', color='black')
        plt.text(len(virus_dict)-0.5, -2, 'Specialist', fontsize=12, ha='left', va='center', color='black')
        plt.grid()
        plt.show()

        # plot host name on x axis and corresponding degree on y axis
        plt.figure(figsize=(10, 6))
        # dot plot for host
        plt.scatter(list(host_dict.keys()), list(host_dict.values()), color='b')
        # turn the x axis labels to be vertical
        plt.xticks(rotation=90)
        plt.title('Host Degree Distribution')
        plt.xlabel('Host')
        plt.ylabel('Degree')
        # Add Generalist and Specialist labels
        plt.text(-0.5, 10, 'Susceptible', fontsize=12, ha='right', va='center', color='black')
        plt.text(len(host_dict)-0.5, 10, 'Resistant', fontsize=12, ha='left', va='center', color='black')
        plt.grid()
        plt.show()

       
    def plot_eigenvector_centrality(self):
        """ Plot the eigenvector centrality of the graph. """
      
        # Virus distribution is green and host distribution is blue
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

        # Virus distribution is green and host distribution is blue
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

        # Virus distribution is green and host distribution is blue
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

    def centrality_boxplot(self):
        """ Plot the centrality measures for host and virus in a boxplot all together. """

        centrality_data = [
            self.eigenvector_virus.values(), self.eigenvector_host.values(),
            self.betweenness_virus.values(), self.betweenness_host.values(),
            self.closeness_virus.values(), self.closeness_host.values()
        ]
        sns.boxplot(data=[list(data) for data in centrality_data],
                palette=["forestgreen", "dodgerblue"] * 3)
        # make a legend where green is virus and blue is host
        plt.legend(["Virus", "Host"], loc='upper left')
        # group the data by the centrality measure
        # first group is eigenvector, second is betweenness, third is closeness
        plt.xticks([0, 2, 4], ['Eigenvector', 'Betweenness', 'Closeness'])
        plt.title('Centrality Measures')

        
    def plot_heatmap(self, prediction_color = "indigo", color_map=["red", "lightpink", "white", "lightskyblue", "blue"], ranges=[0, 0.2, 0.45, 0.55, 0.8, 1]):
        """ Plot the heatmap of the matrix. Let the user choose the colors and ranges of the heatmap.
        Args:
            prediction_color (str): Color for the predictions. Default is "indigo".
            color_map (list): List of colors for the heatmap. Default is ["red", "lightpink", "white", "#a2cffe", "blue"].
            ranges (list): List of ranges for the heatmap. Default is [0, 0.2, 0.45, 0.55, 0.8, 1].
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
        # Add Generalist and Specialist labels
        plt.text(-(len(self.columns)/3), 40, 'Generalist', fontsize=12, color='black')
        plt.text(-(len(self.columns)/3), (len(self.rows)-40), 'Specialist', fontsize=12, ha='left', va='center', color='black')
        # Add Susceptible and Resistant labels
        plt.text(3, (len(self.rows)+(len(self.rows)*0.05)), 'Susceptible', fontsize=12, ha='right', va='center', color='black')
        plt.text(len(self.columns)-2, (len(self.rows)+(len(self.rows)*0.05)), 'Resistant', fontsize=12, ha='left', va='center', color='black')
        # save the figure in the heatmaps folder
        # get matrix name before the first underscore
        plt.savefig('Heatmaps/Heatmap_' + self.title + '_' + 'predictions' if self.probability is False else 'probabiltiies' + '.png')
        plt.show()

    def plot_centrality_time_series(self, measures, centrality_measure:str):
        """ Plot the centrality measures over time. 
        Args:
            measures (list): List of centrality measures for each iteration.
            centrality_measure (str): Name of the centrality measure.
        """
        # Plot a boxplot of the centrality measures for each iteration on one graph
        plt.figure(figsize=(10, 6))
        for i in range(len(measures)):
            plt.boxplot(measures[i], positions=[i])
        plt.title('Centrality Measures Over Time')
        plt.xlabel('Iteration')
        plt.ylabel(centrality_measure)

        
    def plot_mean_centrality_time_series(self, virus_centrality):
        """ Plot the mean centrality measures over time. 
        Args:
            virus_centrality (list): List of centrality measures for each iteration.
        """
        # plot a histogram of the mean centrality measures
        means = []
        for i in range(len(virus_centrality)):
            means.append(np.mean(virus_centrality[i]))
        plt.figure(figsize=(10, 6))
        plt.plot(means)
        plt.title('Mean Centrality Measures Over Time')
        plt.xlabel('Time')
        plt.ylabel('Centrality Measure')
        plt.show()

    def plot_virus_virus_heatmap(self):
        """ Plot the heatmap of the unipartite virus-viruses matrix. """
        # Figure size
        plt.figure(figsize=(10, 10))
        # Mask the diagonal
        mask = np.eye(self.unipartite_viruses.shape[0], dtype=bool)
        # Diagonal is gray and the rest of the matrix is colored using the plasma color map
        cmap = matplotlib.colormaps.get_cmap('plasma')
        cmap.set_bad(color='lightgray')
        sns.heatmap(self.unipartite_viruses, cmap=cmap, mask=mask, cbar_kws={"shrink": .8})
        plt.xlabel("Viruses")
        plt.ylabel("Viruses")
        plt.title(self.title)
        # save the figure in the heatmaps folder
        plt.savefig('Heatmaps/Heatmap_' + self.title + '_unipartite_viruses.png')
        plt.show()
    
    def plot_host_host_heatmap(self):
        """ Plot the heatmap of the unipartite host-host matrix. """  
        # Figure size
        plt.figure(figsize=(10, 10))
        # Mask the diagonal
        mask = np.eye(self.unipartite_hosts.shape[0], dtype=bool)
        # Diagonal is gray and the rest of the matrix is colored using the plasma color map
        cmap = matplotlib.colormaps.get_cmap('plasma')
        cmap.set_bad(color='lightgray')
        sns.heatmap(self.unipartite_hosts, cmap=cmap, mask=mask, cbar_kws={"shrink": .8})
        plt.xlabel("Hosts")
        plt.ylabel("Hosts")
        plt.title(self.title)
        # save the figure in the heatmaps folder
        plt.savefig('Heatmaps/Heatmap_' + self.title + '_unipartite_hosts.png')
        plt.show()
        

    def nestedness_rows(self, pair):
        """Calculate nestedness for rows using the NODF algorithm for the array.
        Args:
            pair (tuple): Tuple containing the indices of the two rows to be compared.
        """
        # generate list of rows to compare
        pair1 = self.input_matrix[pair[0],].astype(int)
        pair2 = self.input_matrix[pair[1],].astype(int)
        return self.compare(pair1, pair2)

    def nestedness_cols(self, pair):
        """Calculate nestedness for cols using the NODF algorithm for the array.
        Args:
            pair (tuple): Tuple containing the indices of the two columns to be compared.
        """
        pair1 = self.input_matrix[:, pair[0]].astype(int)
        pair2 = self.input_matrix[:, pair[1]].astype(int)
        return self.compare(pair1, pair2)

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

    def compare(self, x, y) -> float:
        """Compare two lists containing 0 and 1.
        Args:
            x (list[int]): first list
            y (list[int]): second list
        """
        if sum(x) <= sum(y) or sum(y) == 0:
            return 0
        else:
            counter, total = 0, 0
            # For each pair of values in the two lists, check if we add 1 to the counter and total or just the total
            for i, j in zip(x, y):
                if i == 1 and j == 1:
                    counter += 1
                    total += 1
                elif i == 0 and j == 1:
                    total += 1
        return counter/total * 100
    
    # def compare2(self, x, y) -> float:
    #     """Compare two lists containing 0 and 1.

    #     Args:
    #         x (list[int]): first list
    #         y (list[int]): second list
    #     """

    #     if sum(x) <= sum(y) or sum(y) == 0:
    #         return 0
    #     else:
    #         # the number of -1s in the subtraction contributes to the denominator
    #         count_minus1 = np.count_nonzero(np.subtract(x, y) == -1) 
    #         # the number of 2s in the addition contributes to the numerator and denominator
    #         count_2 = np.count_nonzero(np.add(x, y) == 2) 
    #     return (count_2 / (count_2 + count_minus1)) * 100
    
    def run_nestedness(self):
        """Run the nestedness algorithm for the array. This is not parallelized."""
        # Sort the matrix
        self.sort_matrix()

        # Calculate nestedness for rows
        nrow = []
        for pair in self.pairs(axis=0):
            nrow.append(self.nestedness_rows(pair))
        
        # Calculate nestedness for columns
        ncol = []
        for pair in self.pairs(axis=1):
            ncol.append(self.nestedness_cols(pair))

        nrow = sum(nrow) / (len(self.input_matrix) * (len(self.input_matrix) - 1) / 2)
        ncol = sum(ncol) / (len(self.input_matrix[0]) * (len(self.input_matrix[0]) - 1) / 2)
        return (nrow + ncol) / 2

    def run_parallel(self, num_procs: int = 6):
        """Run multiple process of the compute_feature method.

        This significantly improves run time by using multiple CPU cores.

        Args:
            num_procs (int): Number of cores to be used.
        """

        # The matrix should be sorted for nestedness to work. Will need to re-sort each time?
        self.sort_matrix()

        # Calculate nestedness for rows and columns in parallel
        with Pool() as pool:
            nrow = pool.map(self.nestedness_rows, self.pairs(axis=0))
            ncol = pool.map(self.nestedness_cols, self.pairs(axis=1))

            # Compute the averages
        nrow_avg = sum(nrow) / (len(self.input_matrix) * (len(self.input_matrix) - 1) / 2)
        ncol_avg = sum(ncol) / (len(self.input_matrix[0]) * (len(self.input_matrix[0]) - 1) / 2)
        return (nrow_avg + ncol_avg) / 2
    

    # Calculate modularity
    def calculate_modularity(self):
        """ Calculate the modularity of the graph. """
        self.initialize_graph()
        # Use the Infomap algorithm to calculate modularity
        im = infomap.Infomap()
        im.add_networkx_graph(self.G)
        im.run()
        self.modularity = im.get_modules()
        return self.modularity
    
    # Plot modularity
    def plot_modularity(self, modules):
        """ Plot the modularity of the graph. """
        # Plot the modules
        plt.figure(figsize=(10, 6))
        nx.draw(self.G, node_color=list(modules.values()), with_labels=True)
        plt.title('Modularity')
        plt.show()

    # Connectivity
    def calculate_connectivity(self):
        """ Calculate the connectivity of the graph. """
        self.initialize_graph()
        self.connectivity = nx.average_node_connectivity(self.G)
        return self.connectivity
    
    # Percent edges
    def calculate_percent_edges(self):
        """ Calculate the percent of edges in the graph. """
        # count 1s in the matrix
        return np.count_nonzero(self.input_matrix) / len(self.input_matrix) * len(self.input_matrix[0])
    
    # Clustering coefficient
    def calculate_clustering_coefficient(self):
        """ Calculate the clustering coefficient of the graph. """
        self.initialize_graph()
        self.clustering_coefficient = nx.average_clustering(self.G)
        return self.clustering_coefficient

    # average number of viruses per host
    def calculate_average_viruses_per_host(self):
        """ Calculate the average number of viruses per host. """
        return np.count_nonzero(self.input_matrix) / len(self.input_matrix[0])
    # average number of hosts per virus
    def calculate_average_hosts_per_virus(self):
        """ Calculate the average number of hosts per virus. """
        return np.count_nonzero(self.input_matrix) / len(self.input_matrix)
    