# change directory to one level up
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from VirusHostNetworkAnalysis.prediction_matrix import PredictionMatrix
from VirusHostNetworkAnalysis.null_model import CM

def test_total_degrees():
    test_matrix = PredictionMatrix('tests/test_predictions.tsv')
    test_matrix.make_rectangular_matrix('prediction')
    test_config = CM(test_matrix.virus_host_array)
    degrees = test_config.bootstrap_stats(100)
    assert (sum(degrees[0])) == 8
    assert (sum(degrees[1])) == 8

