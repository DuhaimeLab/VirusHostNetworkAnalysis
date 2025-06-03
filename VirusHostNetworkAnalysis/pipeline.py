from VirusHostNetworkAnalysis.prediction_matrix import PredictionMatrix
from VirusHostNetworkAnalysis.null_model import ConfigurationModel
from VirusHostNetworkAnalysis.null_model import ER
from VirusHostNetworkAnalysis.properties import BipartiteGraph
from tqdm import tqdm
import matplotlib.pyplot as plt
from statistics import mean
import pandas as pd
import os


class Pipeline():
    """ Class to run the pipeline for the VirusHostNetworkAnalysis. 
    
    Args:
    file_path (str): Path to the matrix file.
    null_type (str): Type of null model to use. Options are "ER" or "CM".
    p (float): Probability of edge creation for the ER model.
    num_swaps (int): Number of swaps to perform for the CM model.
    num_runs (int): Number of runs to perform for the null model.
    num_cores (int): Number of cores to use for parallel processing.
    probability (bool): Whether to include the probability matrix or not. Default is False.
    max_iter (int): Maximum number of iterations for the null model. Default is 1000.
    
    """
    def __init__(self, file_path:str, null_type, p:float, cm_method, percent_shuffle_virus:float,
                 percent_shuffle_host:float, num_runs:int, num_cores:int, probability:bool= False, max_iter:int=1000):
        """ Initialize the pipeline with the given parameters. """
        self.file_path = file_path
        self.p = p
        self.num_runs = num_runs
        self.null_type = null_type
        self.num_cores = num_cores
        self.cm_method = cm_method
        self.percent_shuffle_v = percent_shuffle_virus
        self.percent_shuffle_h = percent_shuffle_host

        self.prediction_matrix = PredictionMatrix(self.file_path, False)
        self.prediction_matrix.make_rectangular_matrix()
        if probability is True:
            self.probability_matrix = PredictionMatrix(self.file_path, True)
            self.probability_matrix.make_rectangular_matrix()
            self.probability_properties = BipartiteGraph(self.probability_matrix)

        # Use 10 times the number of edges in the prediction matrix as the number of swaps
        self.num_swaps = self.prediction_matrix.virus_host_array.sum() * 10

    def pipeline_steps_prediction(self):
        """ Run the pipeline for the VirusHostNetworkAnalysis. """
        # PREDICTION MATRIX
        self.prediction_properties = BipartiteGraph(self.prediction_matrix)
        # Calculate the properties
        # degree distribution
        virus_deg, host_deg = self.prediction_properties.calculate_degree()
        self.virus_metrics['degree'].append(virus_deg)
        self.host_metrics['degree'].append(host_deg)

        ### CALCULATE ###
        # Nestedness
        nest = self.prediction_properties.run_parallel(self.num_cores)
        self.nestedness.append(nest)
        print("Nestedness: ", nest)

        # Centrality
        self.run_centrality(self.prediction_properties)

        # Percent edges (# of edges in the matrix / # of possible edges)
        print("Percent edges: ", self.prediction_properties.calculate_percent_edges())

        # Average # of viruses per host
        print("Average number of viruses per host: ", self.prediction_properties.calculate_average_viruses_per_host())

        # Average # of hosts per virus
        print("Average number of hosts per virus: ", self.prediction_properties.calculate_average_hosts_per_virus())

        # Modularity
        modularity = self.prediction_properties.calculate_modularity()
        self.modularity.append(modularity)

    def pipeline_steps_null(self):
        """ Run the pipeline for the null model. Can be ER or CM. """
        if self.null_type == "ER":
            self.null_model = ER(self.prediction_matrix, self.p)
            # Fill the null model
            self.null_model.fill_ER_graph()
        else: #if not ER then will be CM
            self.null_model = ConfigurationModel(self.prediction_matrix)
            # Map methods to their corresponding functions
            cm_methods = {
                "swap": lambda: self.null_model.bootstrap_swaps(self.num_swaps),
                "shuffle": lambda: self.null_model.shuffle_cm(self.percent_shuffle_v, self.percent_shuffle_h),
                "curveball": lambda: self.null_model.curveball_method(self.num_swaps)
            }
            # Execute the selected method
            cm_methods.get(self.cm_method, lambda: None)()
            
        
        # Create the properties object for the null model
        self.null_properties = BipartiteGraph(self.null_model)
        # Centrality
        self.run_centrality(self.null_properties)
        # Nestedness
        nest = self.null_properties.run_parallel(self.num_cores)
        self.nestedness.append(nest)
    
    def run_centrality(self, properties:BipartiteGraph):
        """ Run the centrality calculations for the given properties object. """
        # Calculate centrality
        properties.initialize_graph()
        properties.calculate_centrality(algorithm="eigenvector")
        properties.calculate_centrality(algorithm="betweenness")
        properties.calculate_centrality(algorithm="closeness")

        # Append the centrality values to the metrics lists
        self.virus_metrics["eigenvector"].append(list(properties.eigenvector_virus.values()))
        self.virus_metrics['betweenness'].append(list(properties.betweenness_virus.values()))
        self.virus_metrics['closeness'].append(list(properties.closeness_virus.values()))
        self.host_metrics["eigenvector"].append(list(properties.eigenvector_host.values()))
        self.host_metrics['betweenness'].append(list(properties.betweenness_host.values()))
        self.host_metrics['closeness'].append(list(properties.closeness_host.values()))


    def run_pipeline(self):      
        """ Run the pipeline for the VirusHostNetworkAnalysis. """
        self.virus_metrics = {"eigenvector":[], 
                              "betweenness":[], 
                              "closeness":[],
                              "degree":[]}
        self.host_metrics = {"eigenvector":[],
                             "betweenness":[], 
                             "closeness":[],
                             "degree":[]}
        self.nestedness = []
        self.modularity = []

        # Prediction matrix
        self.pipeline_steps_prediction()

        # Build the unipartite matrix
        self.prediction_properties.unipartite_matrix()

        # If CM is specified, run the configuration model with the given number of swaps
        if self.null_type == "CM" or self.null_type == "ER":
            with tqdm(total=self.num_runs, desc="Running iterations", colour="green") as pbar:
                for i in range(self.num_runs):
                    self.pipeline_steps_null()
                    pbar.update(1)
        
    
    def visualize_scores(self):
        """ Plot the scores of the predictions. """
        self.prediction_matrix.plot_scores()

    def visualize_prediction_heatmap(self, prediction_color:str):
        """ Plot the prediction matrix heatmap. Arguments are passed to the plot_heatmap function.
        
        Args:
            prediction_color (str): Color for the predictions. Default is "indigo".
        """
        self.prediction_properties.plot_heatmap(prediction_color=prediction_color)
    
    def visualize_probability_heatmap(self, color_map, ranges):
        """ Plot the probability matrix heatmap.
        
        Args:
            color_map (str): Color map for the heatmap.
            ranges (list): List of ranges for the color map.
        """
        self.probability_properties.plot_heatmap(color_map=color_map, ranges=ranges)

    def visualize_degree_distribution(self):
        """ Plot the degree distribution of the prediction matrix. """
        self.prediction_properties.plot_degree_distribution()

    def visualize_degree_distribution_null(self):
        """ Plot the degree distribution of the prediction matrix. """
        self.null_properties.plot_degree_distribution()

    def visualize_degree_by_species(self):
        """ Plot the degree by species of the prediction matrix. """
        self.prediction_properties.plot_degree_by_species()
  
    def visualize_prediction_centrality(self):
        """ Plot the centrality of the prediction matrix. """
        # make 3 by 2 figure
        fig, axis = plt.subplots(3, 2, figsize=(14, 16))
        fig.suptitle("Centrality of the Prediction Matrix")
        self.prediction_properties.plot_eigenvector_centrality(axis[0,0], axis[0,1])
        self.prediction_properties.plot_betweenness_centrality(axis[1,0], axis[1,1])
        self.prediction_properties.plot_closeness_centrality(axis[2,0], axis[2,1])

        self.prediction_properties.centrality_boxplot()
        plt.suptitle("Centrality of the Prediction Matrix")


    def visualize_null_centrality(self):
        """ Plot the centrality of the null model. """
        fig, axis = plt.subplots(3, 2, figsize=(14, 16))
        fig.suptitle("Centrality of the Null Model")
        self.null_properties.plot_eigenvector_centrality(axis[0,0], axis[0,1])
        self.null_properties.plot_betweenness_centrality(axis[1,0], axis[1,1])
        self.null_properties.plot_closeness_centrality(axis[2,0], axis[2,1])
        
        self.null_properties.centrality_boxplot()
        plt.suptitle("Centrality of the Null Model")

    def visualize_prediction_vs_null_centrality(self):
        """ Plot the centrality of the prediction matrix vs the null model. """
        self.prediction_properties.plot_prediction_vs_null(self.virus_metrics, self.host_metrics)
    
    def visualize_centrality_over_i(self):
        """ Plot the centrality of the prediction matrix over iterations. 
        Creates different plots for each centrality type and host vs virus. """
        self.null_properties.plot_centrality_time_series(self.virus_metrics['eigenvector'], "Eigenvector Centrality for Virus")
        self.null_properties.plot_centrality_time_series(self.host_metrics['eigenvector'], "Eigenvector Centrality for Host")
        self.null_properties.plot_centrality_time_series(self.virus_metrics['betweenness'], "Betweenness Centrality for Virus")
        self.null_properties.plot_centrality_time_series(self.host_metrics['betweenness'], "Betweenness Centrality for Host")
        self.null_properties.plot_centrality_time_series(self.virus_metrics['closeness'], "Closeness Centrality for Virus")
        self.null_properties.plot_centrality_time_series(self.host_metrics['closeness'], "Closeness Centrality for Host")

    def visualize_unipartite_projection(self):
        """ Plot the unipartite projection of the prediction matrix. """
        self.prediction_properties.unipartite_graph()
        self.prediction_properties.plot_host_host_heatmap()
        self.prediction_properties.plot_virus_virus_heatmap()

    def visualize_nestedness_distribution(self):
        """ Placeholder for nestedness distribution. Need to fix time first. """
        # plot as a distribution
        plt.figure()
        plt.hist(self.nestedness, color='skyblue')
        plt.xlabel("NODF")
        plt.ylabel("Frequency")
        plt.title("NODF Distribution")
        # add a vertical line at the nodf for the prediction matrix
        plt.axvline(x=self.nestedness[0], color='r', linestyle='--', label='Prediction Matrix NODF')

    def export_nestedness(self, directory:str, file_name:str):
        """ Save the nestedness values for each run into a table. 
        
        Args:
            directory (str): Directory to save the file.
            file_name (str): Name of the file to save.
        """
        # Table with nestedness values for each run
        df = pd.DataFrame(columns=["Method", "Iteration", "Percent Shuffle" if self.cm_method == "shuffle" else "Number of Swaps", "Nestedness"])
        df["Method"] = ["Prediction"] + [self.null_type for i in range(self.num_runs)]
        df["Iteration"] = ["Prediction"] + [f"Null iteration {i}" for i in range(1, self.num_runs+1)]
        if self.cm_method == "shuffle":
            df["Percent Shuffle"] = [self.percent_shuffle_v for i in range(self.num_runs+1)]
        else:
            df["Number of Swaps"] = [self.num_swaps for i in range(self.num_runs+1)]
        df["Nestedness"] = self.nestedness

        # Save the dataframe to a csv file in a new folder titled "file_name" in the directory
        df.to_csv(f"{directory}/{file_name}/nestedness.csv", index=False)
    
    def export_centrality(self, directory:str, file_name:str):
        """ Save the centrality values for each run into a table.

        Args:
            directory (str): Directory to save the file.
            file_name (str): Name of the file to save.
        """
        # Table with centrality values. Each run is it's own file and all files are saved in the same directory.
        for i in range(0, self.num_runs+1):
            # Table with centrality values
            df = pd.DataFrame(columns=["Method", "Percent Shuffle" if self.cm_method == "shuffle" else "Number of Swaps", "Node", "Eigenvector Centrality", "Betweenness Centrality", "Closeness Centrality"])
            # Fill in values
            df["Methods"] = ["Prediction"] * (len(self.prediction_matrix.rows) + len(self.prediction_matrix.columns)) if i == 0 else [self.null_type] * (len(self.prediction_matrix.rows) + len(self.prediction_matrix.columns))
            if self.cm_method == "shuffle":
                df["Percent Shuffle"] = [self.percent_shuffle_v for i in range(len(self.prediction_matrix.rows) + len(self.prediction_matrix.columns))]
            else:
                df["Number of Swaps"] = [self.num_swaps for i in range(len(self.prediction_matrix.rows) + len(self.prediction_matrix.columns))]
            # Add the node names to the dataframe
            df["Node"] = (list(self.prediction_matrix.rows) + list(self.prediction_matrix.columns))
            # Add the centrality values to the dataframe
            df["Eigenvector Centrality"] = [val for val in (self.virus_metrics['eigenvector'][i] + self.host_metrics['eigenvector'][i])]
            df["Betweenness Centrality"] = [val for val in (self.virus_metrics['betweenness'][i] + self.host_metrics['betweenness'][i])]
            df["Closeness Centrality"] = [val for val in (self.virus_metrics['closeness'][i] + self.host_metrics['closeness'][i])]
            # Save df to csv
            if i == 0:
                df.to_csv(f"{directory}/{file_name}/centrality_measures_prediction.csv", index=False)
            else:
                # Save the dataframe to a csv file in a new folder titled "file_name" in the directory
                df.to_csv(f"{directory}/{file_name}/centrality_measures_null_{i}.csv", index=False)

    def export_pipeline_data(self, directory:str, file_name:str):
        """ Save the data from the prediction and null models into tables.

        Args:
            directory (str): Directory to save the file.
            file_name (str): Name of the file to save.
        """
        # Create a directory to save the data
        os.makedirs(f"{directory}/{file_name}")

        # Save the nestedness for each run
        self.export_nestedness(directory, file_name)

        # Save the centrality measures for each node and each run
        self.export_centrality(directory, file_name)

        # Make a readme file inside the directory that contains the number of runs, the number of swaps or percent shuffle, the null model type
        with open(f"{directory}/{file_name}/README.txt", "w") as f:
            f.write(f"Number of runs: {self.num_runs}\n")
            f.write(f"Number of swaps: {self.num_swaps}\n")
            f.write(f"Null model type: {self.null_type}\n")
            f.write(f"Percent shuffle: {self.percent_shuffle_v}\n")
            f.write(f"Method: {self.cm_method}\n")
            
            

        
       




