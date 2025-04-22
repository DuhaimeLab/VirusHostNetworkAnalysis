# change directory to one level up
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from VirusHostNetworkAnalysis.prediction_matrix import *
#from VirusHostNetworkAnalysis.prediction_matrix import Calculations

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
               

def test_nestedness():
    test_matrix = [[1, 1, 1, 1, 0], [1, 0, 1, 1, 0], [1, 1, 0, 0, 1], [1, 1, 0, 0, 0], [1, 1, 0, 0, 0]]
    cal = Calculations(test_matrix, True)
    #cal.nestedness()
    assert cal.nestedness() == 58