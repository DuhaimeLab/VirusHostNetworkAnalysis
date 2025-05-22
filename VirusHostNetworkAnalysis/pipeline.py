from VirusHostNetworkAnalysis.prediction_matrix import PredictionMatrix
from VirusHostNetworkAnalysis.null_model import ConfigurationModel
from VirusHostNetworkAnalysis.null_model import ER
from VirusHostNetworkAnalysis.properties import BipartiteGraph
from tqdm import tqdm
import matplotlib.pyplot as plt
from statistics import mean
import pandas as pd


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
    def __init__(self, file_path:str, null_type, p:float, num_swaps:int, num_runs:int, num_cores:int, 
                 probability:bool= False, max_iter:int=1000):
        """ Initialize the pipeline with the given parameters. """
        self.file_path = file_path
        self.p = p
        self.num_swaps = num_swaps
        self.num_runs = num_runs
        self.null_type = null_type
        self.num_cores = num_cores

        self.prediction_matrix = PredictionMatrix(self.file_path, False)
        self.prediction_matrix.make_rectangular_matrix()
        if probability is True:
            self.probability_matrix = PredictionMatrix(self.file_path, True)
            self.probability_matrix.make_rectangular_matrix()
            self.probability_properties = BipartiteGraph(self.probability_matrix)


    def find_num_swaps(self):
        """ Find the number of swaps to use for the null model.
        This is done by running the null model with a number of swaps and plotting the results.
        """
        nodf_swaps = []
        eigen_host = []
        eigen_virus = []
        closeness_host = []
        closeness_virus = []


        prediction_properties = BipartiteGraph(self.prediction_matrix)
        nodf_swaps.append(prediction_properties.run_parallel(self.num_cores))
        cm = ConfigurationModel(self.prediction_matrix)
        for i in range(0, self.num_runs):
            # Perform num_swaps swaps to create the null model
            cm.bootstrap_stats(self.num_swaps)
            # Create the properties object for the null model
            cm_properties = BipartiteGraph(cm)
            # Calculate the properties
            nodf_swaps.append(cm_properties.run_parallel(8))
            cm_properties.calculate_centrality(algorithm="eigenvector")
            closeness_host.append(mean(cm_properties.eigenvector_host.values()))
            closeness_virus.append(mean(cm_properties.eigenvector_virus.values()))
            cm_properties.calculate_centrality(algorithm="closeness")
            eigen_host.append(mean(cm_properties.closeness_host.values()))
            eigen_virus.append(mean(cm_properties.closeness_virus.values()))
            print(i)
        # line graph of the nodf_swaps
        plt.plot(nodf_swaps)
        plt.xlabel(f"Iteration ({self.num_swaps} swaps each)")
        plt.ylabel("NODF")

        # plot as a distribution
        # plt.figure()
        # plt.hist(nodf_swaps, color='skyblue')
        # plt.xlabel("NODF")
        # plt.ylabel("Frequency")
        # plt.title("NODF Distribution")
        # # add a vertical line at the nodf for the prediction matrix
        # plt.axvline(x=nodf_swaps[0], color='r', linestyle='--', label='Prediction Matrix NODF')

        # centrality
        plt.figure()
        plt.plot(eigen_host, label="Host")
        plt.xlabel("Swap #")
        plt.ylabel("Host Eigenvector Centrality")
        # new figure
        plt.figure()
        plt.plot(eigen_virus, label="Virus")
        plt.xlabel("Swap #")
        plt.ylabel("Virus Eigenvector Centrality")

        # plt.figure()
        # plt.hist(eigen_host, color='skyblue')
        # plt.xlabel("Eigenvector Host")
        # plt.ylabel("Frequency")
        # plt.title("Eigenvector Host Distribution")
        # plt.axvline(x=eigen_host[0], color='r', linestyle='--')
        # plt.figure()
        # plt.hist(eigen_virus, color='skyblue')
        # plt.xlabel("Eigenvector Virus")
        # plt.ylabel("Frequency")
        # plt.title("Eigenvector Virus Distribution")
        # plt.axvline(x=eigen_virus[0], color='r', linestyle='--')


        # new figure
        plt.figure()
        plt.plot(closeness_host, label="Host")
        plt.xlabel("Swap #")
        plt.ylabel("Host Closeness Centrality")
        # new figure
        plt.figure()
        plt.plot(closeness_virus, label="Virus")
        plt.xlabel("Swap #")
        plt.ylabel("Virus Closeness Centrality")

        # plt.figure()
        # plt.hist(closeness_host, color='skyblue')
        # plt.xlabel("Closeness Host")
        # plt.ylabel("Frequency")
        # plt.title("Closeness Host Distribution")
        # plt.axvline(x=closeness_host[0], color='r', linestyle='--')
        # plt.figure()
        # plt.hist(closeness_virus, color='skyblue')
        # plt.xlabel("Closeness Virus")
        # plt.ylabel("Frequency")
        # plt.title("Closeness Virus Distribution")
        # plt.axvline(x=closeness_virus[0], color='r', linestyle='--')


        # connectedness

        # modularity

        return

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
        elif self.null_type == "CM":
            self.null_model = ConfigurationModel(self.prediction_matrix)
            # Perform num_swaps swaps to create the null model
            # self.null_model.bootstrap_stats(self.num_swaps)
            self.null_model.method2(0.5, 0.5)
            #self.null_model.curveball_method(self.num_swaps)

        else:
            raise ValueError("Invalid null model type. Choose 'ER' or 'CM'.")
        
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
        if self.null_type == "CM":
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
            color_map (list): List of colors for the heatmap. Default is ["red", "lightpink", "white", "lightskyblue", "blue"].
            ranges (list): List of ranges for the heatmap. Default is [0, 0.2, 0.45, 0.55, 0.8, 1].
        """
        self.prediction_properties.plot_heatmap(prediction_color=prediction_color)
    
    def visualize_probability_heatmap(self, color_map, ranges):
        """ Plot the probability matrix heatmap. """
        self.probability_properties.plot_heatmap(color_map=color_map, ranges=ranges)

    def visualize_degree_distribution(self):
        """ Plot the degree distribution of the prediction matrix. """
        self.prediction_properties.plot_degree_distribution()

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
    
    def visualize_centrality_over_i(self):
        """ Plot the centrality of the prediction matrix over iterations. 
        Creates different plots for each centrality type and host vs virus. """
        self.null_properties.plot_centrality_time_series(self.virus_metrics['eigenvector'], "Eigenvector Centrality")
        self.null_properties.plot_centrality_time_series(self.host_metrics['eigenvector'], "Eigenvector Centrality")
        self.null_properties.plot_centrality_time_series(self.virus_metrics['betweenness'], "Betweenness Centrality")
        self.null_properties.plot_centrality_time_series(self.host_metrics['betweenness'], "Betweenness Centrality")
        self.null_properties.plot_centrality_time_series(self.virus_metrics['closeness'], "Closeness Centrality")
        self.null_properties.plot_centrality_time_series(self.host_metrics['closeness'], "Closeness Centrality")

    def visualize_unipartite_projection(self):
        """ Plot the unipartite projection of the prediction matrix. """
        self.prediction_properties.unipartite_graph()
        self.prediction_properties.plot_host_host_heatmap()
        self.prediction_properties.plot_virus_virus_heatmap()

    def visualize_nestedness_distribution(self):
        """ Placeholder for nestedness distribution. Need to fix time first. """
        return 0

    def export_pipeline_data(self, directory:str, file_name:str):
        """ Save the data from the prediction and null models into tables. """
        # Table with nestedness values
        # Table has nested value in a column
        # Second column first row is "Prediction" and following rows are "Null iteration 1", "Null iteration 2", etc.
        # dataframe object
        df = pd.DataFrame(columns=["Method", "Iteration", "Number of Swaps", "Nestedness"])

        df["Method"] = ["Prediction"] + [self.null_type for i in range(self.num_runs)]
        df["Iteration"] = ["Prediction"] + [f"Null iteration {i}" for i in range(1, self.num_runs+1)]
        df["Number of Swaps"] = [self.num_swaps for i in range(self.num_runs+1)]
        df["Nestedness"] = self.nestedness
        print(df)

        # Save the dataframe to a csv file
        df.to_csv(f"{directory}/{file_name}_nestedness.csv", index=False)

        # Table with centrality values
        df = pd.DataFrame(columns=["Method", "Run", "Number of Swaps", "Node", "Eigenvector Centrality", "Betweenness Centrality", "Closeness Centrality"])
        # Fill in values
        df["Method"] = ["Prediction"] * (len(self.prediction_matrix.rows) + len(self.prediction_matrix.columns)) + \
                   [self.null_type for i in range((self.num_runs) * (len(self.prediction_matrix.rows) + len(self.prediction_matrix.columns)))]
        df["Run"] = ["Prediction"] * (len(self.prediction_matrix.rows) + len(self.prediction_matrix.columns)) + \
                    [f"Null iteration {i}" for i in range(1, self.num_runs+1) for _ in range(len(self.prediction_matrix.rows) + len(self.prediction_matrix.columns))]
        df["Number of Swaps"] = [self.num_swaps for i in range((self.num_runs+1) * (len(self.prediction_matrix.rows) + len(self.prediction_matrix.columns)))]
        # Add the node names to the dataframe
        df["Node"] = (list(self.prediction_matrix.rows) + list(self.prediction_matrix.columns)) * (self.num_runs + 1)
        # Add the centrality values to the dataframe
        df["Eigenvector Centrality"] = [
            val for i in range(self.num_runs+1) 
            for val in (self.virus_metrics['eigenvector'][i] + self.host_metrics['eigenvector'][i])
        ]
        df["Betweenness Centrality"] = [
            val for i in range(self.num_runs+1) 
            for val in (self.virus_metrics['betweenness'][i] + self.host_metrics['betweenness'][i])
        ]
        df["Closeness Centrality"] = [
            val for i in range(self.num_runs+1) 
            for val in (self.virus_metrics['closeness'][i] + self.host_metrics['closeness'][i])
        ]

        # Save df to csv
        df.to_csv(f"{directory}/{file_name}_centrality_measures.csv", index=False)
        # Print first 10 rows of the dataframe
        print(df.head(10))


        
       




