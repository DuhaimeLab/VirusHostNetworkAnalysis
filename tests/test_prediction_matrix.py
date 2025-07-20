# change directory to one level up
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from VirusHostNetworkAnalysis.prediction_matrix import PredictionMatrix
#from VirusHostNetworkAnalysis.prediction_matrix import Calculations

def test_square_matrix():
    """Test that the square matrix is the correct length"""
    test_matrix = PredictionMatrix('tests/test_predictions.tsv')
    test_matrix.get_unique_virus_host()
    assert len(test_matrix.rows) == 6
    assert len(test_matrix.columns) == 3
    # check that there are no duplicates in the unique_viruses and unique_hosts
    assert len(test_matrix.rows) == len(set(test_matrix.rows))
    assert len(test_matrix.columns) == len(set(test_matrix.columns))

def test_initialize_matrix():
    test_matrix = PredictionMatrix('tests/test_predictions.tsv')
    test_matrix.get_unique_virus_host()
    test_matrix.initialize_matrix()
    # check that the matrix is the correct size
    assert len(test_matrix.virus_host_array) == 6
    assert len(test_matrix.virus_host_array[0]) == 3
    # check only filled with 0s
    assert test_matrix.virus_host_array.all() == 0

def test_filling_correctly_pred():
    """Test the matrix is filled correctly"""
    test_matrix = PredictionMatrix('tests/test_predictions.tsv', False)
    # fill matrix for predictions and test that it only has 1s and 0s
    test_matrix.get_unique_virus_host()
    test_matrix.initialize_matrix()
    test_matrix.fill_matrix()
    # assert that matrix contains only 1s and 0s
    assert test_matrix.virus_host_array.all() == 1 or test_matrix.virus_host_array.all() == 0
    # check each spot
    assert test_matrix.virus_host_array[0][0] == 0
    assert test_matrix.virus_host_array[0][1] == 0
    assert test_matrix.virus_host_array[0][2] == 1
    assert test_matrix.virus_host_array[1][0] == 1
    assert test_matrix.virus_host_array[1][1] == 0
    assert test_matrix.virus_host_array[1][2] == 1
    assert test_matrix.virus_host_array[2][0] == 0
    assert test_matrix.virus_host_array[2][1] == 1
    assert test_matrix.virus_host_array[2][2] == 0
    assert test_matrix.virus_host_array[3][0] == 1
    assert test_matrix.virus_host_array[3][1] == 1
    assert test_matrix.virus_host_array[3][2] == 1
    assert test_matrix.virus_host_array[4][0] == 0
    assert test_matrix.virus_host_array[4][1] == 0
    assert test_matrix.virus_host_array[4][2] == 0
    assert test_matrix.virus_host_array[5][0] == 0
    assert test_matrix.virus_host_array[5][1] == 1
    assert test_matrix.virus_host_array[5][2] == 0

def test_filling_correctly_prob():
    """Test the matrix is filled correctly"""
    test_matrix = PredictionMatrix('tests/test_predictions.tsv', True)
    # fill matrix for probabilities and test that it only has floats
    test_matrix.get_unique_virus_host()
    test_matrix.initialize_matrix()
    test_matrix.fill_matrix()
    # assert that matrix contains only floats between 0 and 1
    assert test_matrix.virus_host_array.all() >= 0 and test_matrix.virus_host_array.all() <= 1