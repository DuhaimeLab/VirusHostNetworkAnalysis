# change directory to one level up
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from VirusHostNetworkAnalysis.prediction_matrix import PredictionMatrix
from VirusHostNetworkAnalysis.null_model import ConfigurationModel
from VirusHostNetworkAnalysis.null_model import ER
from VirusHostNetworkAnalysis.properties import BipartiteGraph

def test_ConfigurationModel_total_degrees():
    test_matrix = PredictionMatrix('tests/test_predictions.tsv')
    test_matrix.make_rectangular_matrix()
    test_config = ConfigurationModel(test_matrix)
    test_config.bootstrap_swaps(100)
    test_properties = BipartiteGraph(test_config)
    assert (sum(test_properties.calculate_degree()[0]) == 8)
    assert (sum(test_properties.calculate_degree()[1]) == 8)

def test_ER_total_degrees():
    """Test the total degrees of the ER model."""
    # matrix size is 6x3
    test_matrix = PredictionMatrix('tests/test_predictions.tsv')
    test_matrix.make_rectangular_matrix()
    # 10% probability
    test_er = ER(test_matrix, 0.1)
    test_er.fill_ER_graph()
    print(sum(test_er.virus_host_array))

def test_ConfigurationModel_shuffle_degrees():
    """Test the shuffle degrees of the ConfigurationModel."""
    test_matrix = PredictionMatrix('tests/test_predictions.tsv')
    test_matrix.make_rectangular_matrix()
    test_properties = BipartiteGraph(test_matrix)
    test_config = ConfigurationModel(test_matrix)
    test_config.shuffle_cm(0.8, 0.8)
    test_cm_properties = BipartiteGraph(test_config)
    assert (sorted(test_properties.calculate_degree()) == sorted(test_cm_properties.calculate_degree()))

def test_ConfigurationModel_curveball_degrees():
    """Test the shuffle degrees of the ConfigurationModel."""
    test_matrix = PredictionMatrix('tests/test_predictions.tsv')
    test_matrix.make_rectangular_matrix()
    test_properties = BipartiteGraph(test_matrix)
    test_config = ConfigurationModel(test_matrix)
    test_config.curveball_method(10)
    test_cm_properties = BipartiteGraph(test_config)
    assert (sorted(test_properties.calculate_degree()) == sorted(test_cm_properties.calculate_degree()))
