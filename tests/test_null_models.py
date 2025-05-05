# change directory to one level up
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from VirusHostNetworkAnalysis.prediction_matrix import PredictionMatrix
from VirusHostNetworkAnalysis.null_model import ConfigurationModel
from VirusHostNetworkAnalysis.properties import BipartiteGraph

def test_total_degrees():
    test_matrix = PredictionMatrix('tests/test_predictions.tsv')
    test_matrix.make_rectangular_matrix()
    test_config = ConfigurationModel(test_matrix)
    test_config.bootstrap_stats(100)
    test_properties = BipartiteGraph(test_config)
    assert (test_properties.calculate_degree()[0] == 8)
    assert (test_properties.calculate_degree()[1] == 8)

