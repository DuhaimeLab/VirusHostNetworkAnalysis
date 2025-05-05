# change directory to one level up
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from VirusHostNetworkAnalysis.prediction_matrix import PredictionMatrix
from VirusHostNetworkAnalysis.properties import BipartiteGraph
import numpy as np

def test_sorting():
    """Test the matrix is sorted correctly"""
    test_matrix = PredictionMatrix('tests/test_predictions.tsv', False)
    test_matrix.make_rectangular_matrix()
    test_properties = BipartiteGraph(test_matrix)
    test_properties.sort_matrix()

    # check that sum of each row is the same or smaller than the sum of the row above it
    for i in range(1, len(test_properties.input_matrix)):
        assert sum(test_properties.input_matrix[i]) <= sum(test_properties.input_matrix[i-1])
    # check that sum of each column is the same or smaller than the sum of the column to the left
    for i in range(1, len(test_properties.input_matrix[0])):
        assert sum(test_properties.input_matrix[:, i]) <= sum(test_properties.input_matrix[:, i-1])

    # check that the row and coloumn names are sorted correctly
    assert list(test_properties.rows) == ['v4', 'v2', 'v6', 'v3', 'v1', 'v5'] 
    assert list(test_properties.columns) == ['h3', 'h2', 'h1']

    # check the sums of each row and column
    assert sum(test_properties.input_matrix[0]) == 3
    assert sum(test_properties.input_matrix[1]) == 2
    assert sum(test_properties.input_matrix[2]) == 1
    assert sum(test_properties.input_matrix[3]) == 1
    assert sum(test_properties.input_matrix[4]) == 1
    assert sum(test_properties.input_matrix[5]) == 0
    assert sum(test_properties.input_matrix[:, 0]) == 3
    assert sum(test_properties.input_matrix[:, 1]) == 3
    assert sum(test_properties.input_matrix[:, 2]) == 2


def test_virus_virus():
    """Test the matrix is filled correctly in the host-host portions of the matrix"""
    test_matrix = PredictionMatrix('Sample_Input/Aug4_predictions.tsv')
    test_matrix.make_rectangular_matrix()
    test_properties = BipartiteGraph(test_matrix)
    test_properties.make_square_matrix()
    for i in range(0, len(test_properties.rows)):
        for j in range(len(test_properties.columns), len(test_properties.columns_square)):
            assert test_properties.virus_host_array_square[i][j] == 0
   
def test_host_host():
    """Test the matrix is filled correctly in the virus-virus portions of the matrix"""
    test_matrix = PredictionMatrix('Sample_Input/Aug4_predictions.tsv')
    test_matrix.make_rectangular_matrix()
    test_properties = BipartiteGraph(test_matrix)
    test_properties.make_square_matrix()
    for i in range(len(test_properties.rows), len(test_properties.rows_square)):
        for j in range(0, len(test_properties.columns)):
            assert test_properties.virus_host_array_square[i][j] == 0
                              


# These tests are for nestedness of the first example matrix.
def test_BipartiteGraph_pairs():
    """
    Test the pairs function of the BipartiteGraph class.
    Check the pairs of rows and columns to be compared.
    """

    test_matrix = PredictionMatrix('tests/test_predictions.tsv', False)
    test_matrix.make_rectangular_matrix()
    test_properties = BipartiteGraph(test_matrix)
    test_properties.input_matrix = np.array([[1, 1, 1, 1], 
                                    [1, 1, 1, 0],
                                    [1, 1, 0, 1],
                                    [1, 1, 0, 0],
                                    [1, 1, 0, 0],
                                    [1, 0, 1, 0],
                                    [1, 0, 0, 0]])

    # Test number 1 - 7x4 matrix that appears to be nested
    expected_row_pairs = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
                          (1, 2), (1, 3), (1, 4), (1, 5), (1, 6),
                          (2, 3), (2, 4), (2, 5), (2, 6),
                          (3, 4), (3, 5), (3, 6),
                          (4, 5), (4, 6),
                          (5, 6)]
    expected_col_pairs = [(0, 1), (0, 2), (0, 3), 
                          (1, 2), (1, 3), 
                          (2, 3)]

    assert test_properties.pairs(0) == expected_row_pairs
    assert test_properties.pairs(1) == expected_col_pairs
    assert len(test_properties.pairs(0)) == 21
    assert len(test_properties.pairs(1)) == 6

    # Test number 2 - 7x4 matrix that appears to be modular
    test_properties.input_matrix = np.array([[1, 1, 1, 0], 
                                    [1, 1, 1, 0],
                                    [1, 1, 1, 0],
                                    [1, 1, 1, 0],
                                    [1, 0, 0, 1],
                                    [1, 0, 0, 1],
                                    [1, 0, 0, 1]])
    
    assert test_properties.pairs(0) == expected_row_pairs
    assert test_properties.pairs(1) == expected_col_pairs
    assert len(test_properties.pairs(0)) == 21
    assert len(test_properties.pairs(1)) == 6
    

def test_BipartiteGraph_nestedness_rows():
    """
    Test the nestedness_rows function of the BipartiteGraph class.
    Check each the value for each pair of rows and for the average Nrow.
    """

    # Test number 1 - 7x4 matrix that appears to be nested
    test_matrix = PredictionMatrix('tests/test_predictions.tsv', False)
    test_matrix.make_rectangular_matrix()
    test_properties = BipartiteGraph(test_matrix)
    test_properties.input_matrix = np.array([[1, 1, 1, 1], 
                                    [1, 1, 1, 0],
                                    [1, 1, 0, 1],
                                    [1, 1, 0, 0],
                                    [1, 1, 0, 0],
                                    [1, 0, 1, 0],
                                    [1, 0, 0, 0]])
    expected_val = 78.6
    total = 0
    values = []
    expected_row_pairs = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
                          (1, 2), (1, 3), (1, 4), (1, 5), (1, 6),
                          (2, 3), (2, 4), (2, 5), (2, 6),
                          (3, 4), (3, 5), (3, 6),
                          (4, 5), (4, 6),
                          (5, 6)]
    for expected_row_pair in expected_row_pairs:
        val = test_properties.nestedness_rows(expected_row_pair)
        values.append(val)
        total += val

    assert round(total/21, 1) == expected_val
    assert (values) == [100, 100, 100, 100, 100, 100,
                        0, 100, 100, 100, 100,
                        100, 100, 50, 100,
                        0, 0, 100,
                        0, 100,
                        100]
    
    # Test number 2 - 7x4 matrix that appears to be modular
    test_properties.input_matrix = np.array([[1, 1, 1, 0], 
                                    [1, 1, 1, 0],
                                    [1, 1, 1, 0],
                                    [1, 1, 1, 0],
                                    [1, 0, 0, 1],
                                    [1, 0, 0, 1],
                                    [1, 0, 0, 1]])
    expected_val = 28.6
    total = 0
    values = []
    expected_row_pairs = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
                          (1, 2), (1, 3), (1, 4), (1, 5), (1, 6),
                          (2, 3), (2, 4), (2, 5), (2, 6),
                          (3, 4), (3, 5), (3, 6),
                          (4, 5), (4, 6),
                          (5, 6)]
    for expected_row_pair in expected_row_pairs:
        val = test_properties.nestedness_rows(expected_row_pair)
        values.append(val)
        total += val

    assert round(total/21, 1) == expected_val
    assert (values) == [0, 0, 0, 50, 50, 50,
                        0, 0, 50, 50, 50,
                        0, 50, 50, 50,
                        50, 50, 50,
                        0, 0,
                        0]

def test_BipartiteGraph_nestedness_cols():
    """
    Test the nestedness_cols function of the BipartiteGraph class.
    Check each the value for each pair of columns and for the average Ncol.
    """

    # Test number 1 - 7x4 matrix that appears to be nested
    test_matrix = PredictionMatrix('tests/test_predictions.tsv', False)
    test_matrix.make_rectangular_matrix()
    test_properties = BipartiteGraph(test_matrix)
    test_properties.input_matrix = np.array([[1, 1, 1, 1], 
                                    [1, 1, 1, 0],
                                    [1, 1, 0, 1],
                                    [1, 1, 0, 0],
                                    [1, 1, 0, 0],
                                    [1, 0, 1, 0],
                                    [1, 0, 0, 0]])
    expected_val = 86.1
    total = 0
    values = []
    expected_col_pairs = [(0, 1), (0, 2), (0, 3), 
                          (1, 2), (1, 3), 
                          (2, 3)]
    for expected_col_pair in expected_col_pairs:
        val = test_properties.nestedness_cols(expected_col_pair)
        values.append(round(val, 2))
        total += val

    assert round(total/6, 1) == expected_val
    assert (values) == [100, 100, 100,
                        66.67, 100,
                        50]

    # Test number 2 - 7x4 matrix that appears to be modular
    test_properties.input_matrix = np.array([[1, 1, 1, 0], 
                                    [1, 1, 1, 0],
                                    [1, 1, 1, 0],
                                    [1, 1, 1, 0],
                                    [1, 0, 0, 1],
                                    [1, 0, 0, 1],
                                    [1, 0, 0, 1]])
    expected_val = 50
    total = 0
    values = []
    expected_col_pairs = [(0, 1), (0, 2), (0, 3), 
                          (1, 2), (1, 3), 
                          (2, 3)]
    for expected_col_pair in expected_col_pairs:
        val = test_properties.nestedness_cols(expected_col_pair)
        values.append(round(val, 2))
        total += val

    assert round(total/6, 1) == expected_val
    assert (values) == [100, 100, 100,
                        0, 0,
                        0]

def test_BipartiteGraph_run_parallel():
    """  
    Test the run_parallel function of the BipartiteGraph class. 
    Test that nestedness is calculated correctly for the two test matrices.
    """
    # Test number 1 - 7x4 matrix that appears to be nested
    test_matrix = PredictionMatrix('tests/test_predictions.tsv', False)
    test_matrix.make_rectangular_matrix()
    test_properties = BipartiteGraph(test_matrix)
    test_properties.input_matrix = np.array([[1, 1, 1, 1], 
                                    [1, 1, 1, 0],
                                    [1, 1, 0, 1],
                                    [1, 1, 0, 0],
                                    [1, 1, 0, 0],
                                    [1, 0, 1, 0],
                                    [1, 0, 0, 0]])
    test_properties.rows = np.array(['v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'])
    test_properties.columns = np.array(['h1', 'h2', 'h3', 'h4'])

    expected_val = 82.3
    assert round(test_properties.run_parallel(), 1) == expected_val

    # Test number 2 - 7x4 matrix that appears to be modular
    test_properties.input_matrix = np.array([[1, 1, 1, 0], 
                                    [1, 1, 1, 0],
                                    [1, 1, 1, 0],
                                    [1, 1, 1, 0],
                                    [1, 0, 0, 1],
                                    [1, 0, 0, 1],
                                    [1, 0, 0, 1]])
    
    expected_val = 39.3
    assert round(test_properties.run_parallel(), 1) == expected_val



    