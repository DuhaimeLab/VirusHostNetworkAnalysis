from VirusHostNetworkAnalysis.prediction_matrix import PredictionMatrix
from VirusHostNetworkAnalysis.null_model import ConfigurationModel
from VirusHostNetworkAnalysis.null_model import ER
from VirusHostNetworkAnalysis.properties import BipartiteGraph
from tqdm import tqdm

#GLOBAL VARIABLES
# matrix_file = path
# probability = TRUE or FALSE
# p = 0.5
# num_swaps = 100
# num_runs

class Pipeline():
    """ Class to run the pipeline for the VirusHostNetworkAnalysis. 
    
    Args:
    file_path (str): Path to the matrix file.
    
    
    """
    def __init__(self, file_path:str, null_type, p:float, num_swaps:int, num_runs:int, probability:bool= False, max_iter:int=1000):
        """ Initialize the pipeline with the given parameters. """
        self.file_path = file_path
        self.p = p
        self.num_swaps = num_swaps
        self.null_type = null_type

        self.prediction_matrix = PredictionMatrix(self.file_path, probability)
        self.prediction_matrix.make_rectangular_matrix()

    def pipeline_steps(self):
        """ Run the pipeline for the VirusHostNetworkAnalysis. """
        # Create the null model
        # if self.null_type == "ER":
        #     self.null_model = ER(self.prediction_matrix, self.p)
        #     # Fill the null model
        #     self.null_model.fill_ER_graph()
        # elif self.null_type == "CM":
        self.null_model = ConfigurationModel(self.prediction_matrix)
            # Fill the null model
        self.null_model.bootstrap_stats(self.num_swaps)
        # elif self.null_type == "None":
        #     pass
        # else:
        #     raise ValueError("Invalid null model type. Choose 'ER' or 'CM'.")
    
        
        # Create the properties object
        if self.null_type == "None":
            self.properties = BipartiteGraph(self.prediction_matrix)
        else:
            self.properties = BipartiteGraph(self.null_model)
        
        # Calculate the properties
        # Degree distribution
        virus_deg, host_deg = self.properties.calculate_degree()
        #self.properties.plot_degree_distribution()
        self.virus_metrics['degree'].append(virus_deg)
        self.host_metrics['degree'].append(host_deg)


        # Centrality
        self.properties.initialize_graph()
        self.properties.calculate_centrality(algorithm="eigenvector")
        #self.properties.calculate_centrality(algorithm="betweenness")
        #self.properties.calculate_centrality(algorithm="closeness")
        self.virus_metrics["eigenvector"].append(list(self.properties.eigenvector_virus.values()))
        # self.virus_metrics['betweenness'].append(list(self.properties.betweenness_virus.values()))
        # self.virus_metrics['closeness'].append(list(self.properties.closeness_virus.values()))
        self.host_metrics["eigenvector"].append(list(self.properties.eigenvector_host.values()))
        # self.host_metrics['betweenness'].append(list(self.properties.betweenness_host.values()))
        # self.host_metrics['closeness'].append(list(self.properties.closeness_host.values()))



    def run_pipeline(self, num_runs:int):      
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
        self.properties = BipartiteGraph(self.prediction_matrix)

        # Centrality
        self.properties.initialize_graph()
        self.properties.calculate_centrality(algorithm="eigenvector")
        #self.properties.calculate_centrality(algorithm="betweenness")
        #self.properties.calculate_centrality(algorithm="closeness")
        self.virus_metrics["eigenvector"].append(list(self.properties.eigenvector_virus.values()))
        # self.virus_metrics['betweenness'].append(list(self.properties.betweenness_virus.values()))
        # self.virus_metrics['closeness'].append(list(self.properties.closeness_virus.values()))
        self.host_metrics["eigenvector"].append(list(self.properties.eigenvector_host.values()))
        # self.host_metrics['betweenness'].append(list(self.properties.betweenness_host.values()))
        # self.host_metrics['closeness'].append(list(self.properties.closeness_host.values()))

        # if self.null_type == "None":
        #     self.properties.plot_centrality_time_series(self.virus_metrics["eigenvector"])
        #     self.properties.plot_mean_centrality_time_series(self.virus_metrics["eigenvector"])
        #     return
        if self.null_type == "CM":
            with tqdm(total=num_runs, desc="Running iterations", colour="green") as pbar:
                for i in range(num_runs):
                    self.pipeline_steps()
                    pbar.update(1)

            self.properties.plot_centrality_time_series(self.virus_metrics["eigenvector"])
            self.properties.plot_mean_centrality_time_series(self.virus_metrics["eigenvector"])