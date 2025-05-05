from VirusHostNetworkAnalysis.prediction_matrix import PredictionMatrix
from VirusHostNetworkAnalysis.null_model import ConfigurationModel
from VirusHostNetworkAnalysis.null_model import ER
from VirusHostNetworkAnalysis.properties import BipartiteGraph

#GLOBAL VARIABLES
# matrix_file = path
# probability = TRUE or FALSE
# p = 0.5
# num_iterations = 100

class Pipeline():
    """ Class to run the pipeline for the VirusHostNetworkAnalysis. 
    
    Args:
    file_path (str): Path to the matrix file.
    
    
    """
    def __init__(self, file_path:str, model_type, p:float, num_iterations:int, probability:bool= False):
        self.file_path = file_path
        self.p = p
        self.num_iterations = num_iterations

        self.prediction_matrix = PredictionMatrix(self.file_path, probability)
        self.prediction_matrix.make_rectangular_matrix()

