# change directory to one level up
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from VirusHostNetworkAnalysis.prediction_matrix import PredictionMatrix


def test_square_matrix():
    """Test that the square matrix is the correct length"""
    test_matrix = PredictionMatrix('tests/test_predictions.tsv')
    test_matrix.get_unique_virus_host()
    assert len(test_matrix.unique_viruses) == 2
    assert len(test_matrix.unique_hosts) == 17

def test_filling_correctly_pred():
    """Test the matrix is filled correctly"""
    test_matrix = PredictionMatrix('tests/test_predictions.tsv')
    # fill matrix for predictions and test that it only has 1s and 0s
    test_matrix.get_unique_virus_host()
    test_matrix.initialize_matrix('prediction')
    test_matrix.fill_matrix('prediction')
    # assert that matrix contains only 1s and 0s
    assert test_matrix.virus_host_array.all() == 1 or test_matrix.virus_host_array.all() == 0

def test_filling_correctly_prob():
    """Test the matrix is filled correctly"""
    test_matrix = PredictionMatrix('tests/test_predictions.tsv')
    # fill matrix for probabilities and test that it only has floats
    test_matrix.get_unique_virus_host()
    test_matrix.initialize_matrix('probability')
    test_matrix.fill_matrix('probability')
    # assert that matrix contains only floats between 0 and 1
    assert test_matrix.virus_host_array.all() >= 0 and test_matrix.virus_host_array.all() <= 1

def test_sorting():
    """Test the matrix is sorted correctly"""
    test_matrix = PredictionMatrix('tests/test_predictions.tsv')
    test_matrix.get_unique_virus_host()
    test_matrix.initialize_matrix('prediction')
    test_matrix.fill_matrix('prediction')
    test_matrix.sort_matrix()
    # check that the first element is 1 and the last element is 0
    assert test_matrix.virus_host_array[0][0] == 1
    assert test_matrix.virus_host_array[-1][-1] == 0


