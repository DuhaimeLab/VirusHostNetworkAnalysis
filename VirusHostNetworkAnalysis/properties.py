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
import itertools


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
    
    def sort_rows_cols(self, axis:int, bipartite:bool = True):
        """ Find the sum of each row and move rows with the highest sum to the top of the matrix. 
        
        Args:
            axis (int): Axis to be used when sorting the matrix. 0 for rows, 1 for columns.
        
        """
        if bipartite is True:
            # Calculate the sum of each row or column
            counts = np.sum(self.input_matrix, axis=axis)
            # Sort inidices based on the counts. This method sorts the sums in descending order, so any ties will be broken by descending order of the index.
            sorted_indices = np.argsort(counts)[::-1]
            # Sort the matrix and row names based on the sorted indices
            if axis == 1:
                self.rows = self.rows[sorted_indices]
                self.input_matrix = self.input_matrix[sorted_indices]
            elif axis == 0:
                self.columns = self.columns[sorted_indices]
                self.input_matrix = self.input_matrix[:, sorted_indices]
        else:
            # Calculate the sum of each row or column
            counts_virus = np.sum(self.unipartite_viruses, axis=axis)
            counts_host = np.sum(self.unipartite_hosts, axis=axis)
            # Sort inidices based on the counts. This method sorts the sums in descending order, so any ties will be broken by descending order of the index.
            sorted_indices_virus = np.argsort(counts_virus)[::-1]
            sorted_indices_host = np.argsort(counts_host)[::-1]
            # Sort the matrices based on the sorted indices
            if axis == 1:
                self.unipartite_viruses = self.unipartite_viruses[sorted_indices_virus]
                self.unipartite_hosts = self.unipartite_hosts[sorted_indices_host]
            elif axis == 0:
                self.unipartite_viruses = self.unipartite_viruses[:, sorted_indices_virus]
                self.unipartite_hosts = self.unipartite_hosts[:, sorted_indices_host]


    def sort_matrix(self, bipartite:bool = True):
        """ Sort the matrix by rows and columns. """
        # Sort the rows
        self.sort_rows_cols(1, bipartite)
        # Sort the columns
        self.sort_rows_cols(0, bipartite)


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
    
    def make_square_matrix(self):
        """ Expand the matrix to make it square by adding rows for the hosts and columns for the viruses. """
        # rectangular matrix should be sorted first
        self.sort_matrix(True)
        # add number of hosts to the rows
        self.virus_host_array_square = np.concatenate((self.input_matrix, np.zeros((len(self.columns), len(self.columns)), dtype=bool if self.probability is False else float)), axis=0)
        self.rows_square = np.concatenate((self.rows, self.columns))
        # add number of viruses to the columns
        self.virus_host_array_square = np.concatenate((self.virus_host_array_square, np.zeros((len(self.rows_square), len(self.rows)), bool if self.probability is False else float)), axis=1)
        self.columns_square = np.concatenate((self.columns, self.rows))
        # Fill with 1s
        self.fill_bottom_right()

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
        
        # Sort the unipartite matrices
        self.sort_matrix(False)
        
        

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
                
    def chunks(self, l, n):
        """Divide a list of nodes `l` in `n` chunks"""
        l_c = iter(l)
        while 1:
            x = tuple(itertools.islice(l_c, n))
            if not x:
                return
            yield x

    # Calculate the centrality of the graph
    def calculate_centrality(self, algorithm = "eigenvector", max_iter = 1000):
        """ Calculate the centrality of the graph. """
        self.initialize_graph()
        # lower the algorithm name to make it case insensitive
        algorithm = algorithm.lower()
        if algorithm == "eigenvector":
            self.eigenvector = nx.eigenvector_centrality(self.G, max_iter=max_iter)
            # Split the eigenvector centrality into virus and host
            self.eigenvector_virus = {k: self.eigenvector[k] for k in self.rows if k in self.eigenvector}
            self.eigenvector_host = {k: self.eigenvector[k] for k in self.columns if k in self.eigenvector}
        elif algorithm == "betweenness":
            # self.betweenness = nx.betweenness_centrality(self.G)
            # # Split the betweenness centrality into virus and host
            # self.betweenness_virus = {k: self.betweenness[k] for k in self.rows if k in self.betweenness}
            # self.betweenness_host = {k: self.betweenness[k] for k in self.columns if k in self.betweenness}
    
            """Parallel betweenness centrality function"""
            p = Pool(processes=None)
            from multiprocessing import cpu_count
            node_divisor = cpu_count() * 4
            node_chunks = list(self.chunks(self.G.nodes(), 1 if (self.G.order() // node_divisor)==0 else self.G.order() // node_divisor))
            num_chunks = len(node_chunks)
            bt_sc = p.starmap(
                nx.betweenness_centrality_subset,
                zip(
                    [self.G] * num_chunks,
                    node_chunks,
                    [list(self.G)] * num_chunks,
                    [True] * num_chunks,
                    [None] * num_chunks,
                ),
            )

            # Reduce the partial solutions
            from collections import defaultdict
            bt_c = defaultdict(float, bt_sc[0])
            for bt in bt_sc[1:]:
                for n in bt:
                    bt_c[n] += bt[n]
            
            # split the betweenness centrality into virus and host
            self.betweenness_virus = {k: bt_c[k] for k in self.rows if k in bt_c}
            self.betweenness_host = {k: bt_c[k] for k in self.columns if k in bt_c}
            

            
        elif algorithm == "closeness":
            self.closeness = nx.closeness_centrality(self.G)
            # Split the closeness centrality into virus and host
            self.closeness_virus = {k: self.closeness[k] for k in self.rows if k in self.closeness}
            self.closeness_host = {k: self.closeness[k] for k in self.columns if k in self.closeness}   
        else:
            raise ValueError("Algorithm not supported. Choose from 'eigenvector', 'betweenness', or 'closeness'.")

    def calculate_degree(self):
        """ Calculate the degree of each node in the graph. 
        
        Returns:
            list: A list of two lists containing degrees, one for virus and one for host.
        """
        
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

        sns.kdeplot(host_degree, color='dodgerblue', fill=True)
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

        # put next to each other
        fig, axis = plt.subplots(1, 2, figsize=(20, 6), gridspec_kw={'wspace': 0.3})
        # dot plot for virus
        axis[0].plot(list(virus_dict.keys()), list(virus_dict.values()), color='g')
        # remove x axis labels
        axis[0].set_xticks([])
        axis[0].set_title('Virus Degree Distribution')
        axis[0].set_xlabel('Virus')
        axis[0].set_ylabel('Degree')
        # Add Generalist and Specialist labels
        axis[0].text(-0.5, -2, 'Generalist', fontsize=12, ha='right', va='center', color='black')
        axis[0].text(len(virus_dict)-0.5, -2, 'Specialist', fontsize=12, ha='left', va='center', color='black')
        axis[0].grid()

        # plot host name on x axis and corresponding degree on y axis
        # dot plot for host
        axis[1].plot(list(host_dict.keys()), list(host_dict.values()), color='dodgerblue', marker='o', linestyle='-')
        # turn the x axis labels to be vertical
        axis[1].set_xticks([])
        axis[1].set_title('Host Degree Distribution')
        axis[1].set_xlabel('Host')
        axis[1].set_ylabel('Degree')
        # Add Generalist and Specialist labels
        axis[1].text(-0.5, 0, 'Susceptible', fontsize=12, ha='right', va='center', color='black')
        axis[1].text(len(host_dict)-0.5, -10, 'Resistant', fontsize=12, ha='left', va='center', color='black')
        axis[1].grid()

       
    def plot_eigenvector_centrality(self, ax1, ax2):
        """ Plot the eigenvector centrality of the graph. """
      
        # Virus distribution is green and host distribution is blue
        # Plot for the virus eigenvector centrality
        sns.kdeplot(list(self.eigenvector_virus.values()), color='g', fill=True, ax=ax1)
        ax1.ticklabel_format(style='plain')
        ax1.set_title('Eigenvector Centrality for Viruses')
        ax1.set_xlabel('Eigenvector Centrality')
        ax1.set_ylabel('Frequency')
        ax1.grid()

        # Plot for the host eigenvector centrality
        sns.kdeplot(list(self.eigenvector_host.values()), color='dodgerblue', fill=True, ax=ax2)
        ax2.ticklabel_format(style='plain')
        ax2.set_title('Eigenvector Centrality for Hosts')
        ax2.set_xlabel('Eigenvector Centrality')
        ax2.set_ylabel('Frequency')
        ax2.grid()

    def plot_betweenness_centrality(self, ax1, ax2):
        """ Plot the betweenness centrality of the graph. """

        # Virus distribution is green and host distribution is blue
        # Plot for the virus betweenness centrality
        sns.kdeplot(list(self.betweenness_virus.values()), color='g', fill=True, ax=ax1)
        ax1.ticklabel_format(style='plain')
        ax1.set_title('Betweenness Centrality for Viruses')
        ax1.set_xlabel('Betweenness Centrality')
        ax1.set_ylabel('Frequency')
        ax1.grid()

        # Plot for the host betweenness centrality
        sns.kdeplot(list(self.betweenness_host.values()), color='dodgerblue', fill=True, ax=ax2)
        ax2.ticklabel_format(style='plain')
        ax2.set_title('Betweenness Centrality for Hosts')
        ax2.set_xlabel('Betweenness Centrality')
        ax2.set_ylabel('Frequency')
        ax2.grid()

    def plot_closeness_centrality(self, ax1, ax2):
        """ Plot the closeness centrality of the graph. """

        # Virus distribution is green and host distribution is blue
        # Plot for the virus closeness centrality
        sns.kdeplot(list(self.closeness_virus.values()), color='g', fill=True, ax=ax1)
        ax1.ticklabel_format(style='plain')
        ax1.set_title('Closeness Centrality for Viruses')
        ax1.set_xlabel('Closeness Centrality')
        ax1.set_ylabel('Frequency')
        ax1.grid()

        # Plot for the host closeness centrality
        sns.kdeplot(list(self.closeness_host.values()), color='dodgerblue', fill=True, ax=ax2)
        ax2.ticklabel_format(style='plain')
        ax2.set_title('Closeness Centrality for Hosts')
        ax2.set_xlabel('Closeness Centrality')
        ax2.set_ylabel('Frequency')
        ax2.grid()

    def centrality_boxplot(self):
        """ Plot the centrality measures for host and virus in a boxplot all together. """
        fig, axis = plt.subplots(1, 3, figsize=(15, 6))
        sns.boxplot(data=[list(data) for data in [self.eigenvector_virus.values(), self.eigenvector_host.values()]],
                palette=["forestgreen", "dodgerblue"], ax=axis[0])
        sns.boxplot(data=[list(data) for data in [self.betweenness_virus.values(), self.betweenness_host.values()]],
                palette=["forestgreen", "dodgerblue"], ax=axis[1])
        sns.boxplot(data=[list(data) for data in [self.closeness_virus.values(), self.closeness_host.values()]],
            palette=["forestgreen", "dodgerblue"], ax=axis[2])
        axis[0].set_ylabel("Centrality")
        axis[0].set_title("Eigenvector Centrality")
        axis[1].set_title("Betweenness Centrality")
        axis[2].set_title("Closeness Centrality")
        axis[0].set_xticks([0, 1])
        axis[0].set_xticklabels(["Virus", "Host"])
        axis[1].set_xticks([0, 1])
        axis[1].set_xticklabels(["Virus", "Host"])
        axis[2].set_xticks([0, 1])
        axis[2].set_xticklabels(["Virus", "Host"])

    def plot_prediction_vs_null(self, virus_metrics, host_metrics):
        """ Plot the centrality measurements for the predictions and the null model. 
        
        Args:
            virus_metrics (dict): Dictionary of centrality measures for viruses. Includes eigenvector, betweenness, and closeness.
            host_metrics (dict): Dictionary of centrality measures for hosts. Includes eigenvector, betweenness, and closeness.
        """
        # boxplot for the viruses
        # 3 images in one figure, one for each centrality measure
        fig, axis = plt.subplots(1, 3, figsize=(15, 6))
        # Set title to be "Centrality Measures for Viruses"
        fig.suptitle("Centrality Measures for Viruses", fontsize=16)
        sns.boxplot(data=[list(data) for data in [virus_metrics["eigenvector"][0],
                                                  virus_metrics["eigenvector"][-1]]], color="forestgreen", ax=axis[0])
        axis[0].set_title("Eigenvector Centrality")
        # set labels 
        axis[0].set_xticks([0, 1])
        axis[0].set_xticklabels(["Predicted", "Null"])
        sns.boxplot(data=[list(data) for data in [virus_metrics["betweenness"][0],
                                                  virus_metrics["betweenness"][-1]]], color="forestgreen", ax=axis[1])
        axis[1].set_title("Betweenness Centrality")
        axis[1].set_xticks([0, 1])
        axis[1].set_xticklabels(["Predicted", "Null"])
        sns.boxplot(data=[list(data) for data in [virus_metrics["closeness"][0],
                                                  virus_metrics["closeness"][-1]]], color="forestgreen", ax=axis[2])
        axis[2].set_title("Closeness Centrality")
        axis[2].set_xticks([0, 1])
        axis[2].set_xticklabels(["Predicted", "Null"])

        # boxplot for the hosts
        # 3 images in one figure, one for each centrality measure
        fig, axis = plt.subplots(1, 3, figsize=(15, 6))
        # Set title to be "Centrality Measures for Hosts"
        fig.suptitle("Centrality Measures for Hosts", fontsize=16)
        sns.boxplot(data=[list(data) for data in [host_metrics["eigenvector"][0],
                                                  host_metrics["eigenvector"][-1]]], color="dodgerblue", ax=axis[0])  
        axis[0].set_title("Eigenvector Centrality")
        axis[0].set_xticks([0, 1])
        axis[0].set_xticklabels(["Predicted", "Null"])
        sns.boxplot(data=[list(data) for data in [host_metrics["betweenness"][0],
                                                  host_metrics["betweenness"][-1]]], color="dodgerblue", ax=axis[1])
        axis[1].set_title("Betweenness Centrality")
        axis[1].set_xticks([0, 1])
        axis[1].set_xticklabels(["Predicted", "Null"])
        sns.boxplot(data=[list(data) for data in [host_metrics["closeness"][0],
                                                  host_metrics["closeness"][-1]]], color="dodgerblue", ax=axis[2])
        axis[2].set_title("Closeness Centrality")
        axis[2].set_xticks([0, 1])
        axis[2].set_xticklabels(["Predicted", "Null"])


    def plot_heatmap(self, prediction_color = "indigo", color_map=["red", "lightpink", "white", "lightskyblue", "blue"], ranges=[0, 0.2, 0.45, 0.55, 0.8, 1]):
        """ Plot the heatmap of the matrix. Let the user choose the colors and ranges of the heatmap.
        Args:
            prediction_color (str): Color for the predictions. Default is "indigo".
            color_map (list): List of colors for the heatmap. Default is ["red", "lightpink", "white", "#a2cffe", "blue"].
            ranges (list): List of ranges for the heatmap. Default is [0, 0.2, 0.45, 0.55, 0.8, 1].
        """

        # Sort the matrix so that patterns are visible
        self.sort_matrix(True)

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
        # title variable is the name of the predictions file without the .tsv
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
            # fill the boxplot with a color
            plt.boxplot(measures[i], positions=[i], patch_artist=True, boxprops=dict(facecolor='aliceblue', color='black'), medianprops=dict(color='black'))
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
        # Sort the matrix so that patterns are visible
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
        Returns:
            float: Nestedness value for the pair of rows.
        """
        # generate list of rows to compare
        pair1 = self.input_matrix[pair[0],].astype(int)
        pair2 = self.input_matrix[pair[1],].astype(int)
        return self.compare(pair1, pair2)

    def nestedness_cols(self, pair):
        """Calculate nestedness for cols using the NODF algorithm for the array.
        Args:
            pair (tuple): Tuple containing the indices of the two columns to be compared.
        Returns:
            float: Nestedness value for the pair of columns.
        """
        pair1 = self.input_matrix[:, pair[0]].astype(int)
        pair2 = self.input_matrix[:, pair[1]].astype(int)
        return self.compare(pair1, pair2)

    def pairs(self, axis: int = 0) -> List[Tuple[int, int]]:
        """Determine all possible i-j pairs.

        Args:
            axis (int): Axis to be used when determining all pairs.
        Returns:
            List[Tuple[int, int]]: List of tuples containing all possible i-j pairs. 
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
        Returns:
            float: Nestedness value for the pair of lists.
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
    
    def run_nestedness(self):
        """Run the nestedness algorithm for the array. This is not parallelized.
        
        Returns:
            float: Nestedness value for the array.
        """
        # Sort the matrix
        self.sort_matrix(True)

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
        Returns:
            float: Nestedness value for the array.
        """

        # The matrix should be sorted for nestedness to work. Will need to re-sort each time
        self.sort_matrix(True)

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
        """ Calculate the modularity of the graph. 
        
        Returns:
            modularity (float): The modularity of the graph.
        """
        self.initialize_graph()
        # Use the Infomap algorithm to calculate modularity
        im = infomap.Infomap()
        im.add_networkx_graph(self.G)
        im.run()
        self.modularity = im.get_modules()
        # get average modularity
        self.modularity = np.mean(list(self.modularity.values()))
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
        """ Calculate the connectivity of the graph. 
        
        Returns:
            connectivity (float): The connectivity of the graph.
        """
        self.initialize_graph()
        self.connectivity = nx.average_node_connectivity(self.G)
        return self.connectivity
    
    # Percent edges
    def calculate_percent_edges(self):
        """ Calculate the percent of edges in the graph. 

        Returns:
            float: The percent of edges in the graph.
        """
        # count 1s in the matrix
        return np.count_nonzero(self.input_matrix) / (len(self.input_matrix) * len(self.input_matrix[0]))
    
    # Clustering coefficient
    def calculate_clustering_coefficient(self):
        """ Calculate the clustering coefficient of the graph. 
        
        Return:
            clustering_coefficient (float): The clustering coefficient of the graph.
        """
        self.initialize_graph() #initialize the graph G
        self.clustering_coefficient = nx.average_clustering(self.G)
        return self.clustering_coefficient

    # average number of viruses per host
    def calculate_average_viruses_per_host(self):
        """ Calculate the average number of viruses per host. 
        
        Returns:
            float: The average number of viruses per host.
        """
        return np.count_nonzero(self.input_matrix) / len(self.input_matrix[0])
    # average number of hosts per virus
    def calculate_average_hosts_per_virus(self):
        """ Calculate the average number of hosts per virus. 
        
        Returns:
            float: The average number of hosts per virus.
        """
        return np.count_nonzero(self.input_matrix) / len(self.input_matrix)
    