# change directory to one level up
import os
import sys

import test
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from VirusHostNetworkAnalysis.prediction_matrix import PredictionMatrix
#from VirusHostNetworkAnalysis.prediction_matrix import Calculations

def test_square_matrix():
    """Test that the square matrix is the correct length"""
    test_matrix = PredictionMatrix('tests/test_predictions.tsv')
    test_matrix.get_unique_virus_host()
    assert len(test_matrix.unique_viruses) == 6
    assert len(test_matrix.unique_hosts) == 3
    # check that there are no duplicates in the unique_viruses and unique_hosts
    assert len(test_matrix.unique_viruses) == len(set(test_matrix.unique_viruses))
    assert len(test_matrix.unique_hosts) == len(set(test_matrix.unique_hosts))

def test_initialize_matrix():
    test_matrix = PredictionMatrix('tests/test_predictions.tsv')
    test_matrix.get_unique_virus_host()
    test_matrix.initialize_matrix('prediction')
    # check that the matrix is the correct size
    assert len(test_matrix.virus_host_array) == 6
    assert len(test_matrix.virus_host_array[0]) == 3
    # check only filled with 0s
    assert test_matrix.virus_host_array.all() == 0

def test_filling_correctly_pred():
    """Test the matrix is filled correctly"""
    test_matrix = PredictionMatrix('tests/test_predictions.tsv')
    # fill matrix for predictions and test that it only has 1s and 0s
    test_matrix.get_unique_virus_host()
    test_matrix.initialize_matrix('prediction')
    test_matrix.fill_matrix('prediction')
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

    # check that sum of each row is the same or smaller than the sum of the row above it
    for i in range(1, len(test_matrix.virus_host_array)):
        assert sum(test_matrix.virus_host_array[i]) <= sum(test_matrix.virus_host_array[i-1])
    # check that sum of each column is the same or smaller than the sum of the column to the left
    for i in range(1, len(test_matrix.virus_host_array[0])):
        assert sum(test_matrix.virus_host_array[:, i]) <= sum(test_matrix.virus_host_array[:, i-1])

    # check that the row and coloumn names are sorted correctly
    assert list(test_matrix.rows) == ['v4', 'v2', 'v6', 'v3', 'v1', 'v5'] 
    assert list(test_matrix.columns) == ['h3', 'h2', 'h1']

    # check the sums of each row and column
    assert sum(test_matrix.virus_host_array[0]) == 3
    assert sum(test_matrix.virus_host_array[1]) == 2
    assert sum(test_matrix.virus_host_array[2]) == 1
    assert sum(test_matrix.virus_host_array[3]) == 1
    assert sum(test_matrix.virus_host_array[4]) == 1
    assert sum(test_matrix.virus_host_array[5]) == 0
    assert sum(test_matrix.virus_host_array[:, 0]) == 3
    assert sum(test_matrix.virus_host_array[:, 1]) == 3
    assert sum(test_matrix.virus_host_array[:, 2]) == 2


def test_host_host():
    """Test the matrix is filled correctly in the host-host portions of the matrix"""
    test_matrix = PredictionMatrix('Sample_Input/Aug4_predictions.tsv')
    matrix_square = test_matrix.make_square_matrix('prediction')
    for i in range(0, len(test_matrix.unique_viruses)):
        for j in range(len(test_matrix.unique_hosts), len(test_matrix.columns_square)):
            assert matrix_square[i][j] == 0
   
def test_virus_virus():
    """Test the matrix is filled correctly in the virus-virus portions of the matrix"""
    test_matrix = PredictionMatrix('Sample_Input/Aug4_predictions.tsv')
    matrix_square = test_matrix.make_square_matrix('prediction')
    for i in range(len(test_matrix.unique_viruses), len(test_matrix.rows_square)):
        for j in range(0, len(test_matrix.unique_hosts)):
            assert matrix_square[i][j] == 0
                              
