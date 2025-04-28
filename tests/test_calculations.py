# change directory to one level up
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from VirusHostNetworkAnalysis.prediction_matrix import Calculations


# These tests are for nestedness of the first example matrix.
def test_Calculations_pairs():
    """
    Test the pairs function of the Calculations class.
    Check the pairs of rows and columns to be compared.
    """

    # Test number 1 - 7x4 matrix that appears to be nested
    test_matrix = [[1, 1, 1, 1], 
                   [1, 1, 1, 0],
                   [1, 1, 0, 1],
                   [1, 1, 0, 0],
                   [1, 1, 0, 0],
                   [1, 0, 1, 0],
                   [1, 0, 0, 0]]
    expected_row_pairs = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
                          (1, 2), (1, 3), (1, 4), (1, 5), (1, 6),
                          (2, 3), (2, 4), (2, 5), (2, 6),
                          (3, 4), (3, 5), (3, 6),
                          (4, 5), (4, 6),
                          (5, 6)]
    expected_col_pairs = [(0, 1), (0, 2), (0, 3), 
                          (1, 2), (1, 3), 
                          (2, 3)]

    assert Calculations(test_matrix, True).pairs(0) == expected_row_pairs
    assert Calculations(test_matrix, True).pairs(1) == expected_col_pairs
    assert len(Calculations(test_matrix, True).pairs(0)) == 21
    assert len(Calculations(test_matrix, True).pairs(1)) == 6

    # Test number 2 - 7x4 matrix that appears to be modular
    test_matrix_2 = [[1, 1, 1, 0], 
                   [1, 1, 1, 0],
                   [1, 1, 1, 0],
                   [1, 1, 1, 0],
                   [1, 0, 0, 1],
                   [1, 0, 0, 1],
                   [1, 0, 0, 1]]
    
    assert Calculations(test_matrix_2, True).pairs(0) == expected_row_pairs
    assert Calculations(test_matrix_2, True).pairs(1) == expected_col_pairs
    assert len(Calculations(test_matrix_2, True).pairs(0)) == 21
    assert len(Calculations(test_matrix_2, True).pairs(1)) == 6
    

def test_Calculations_nestedness_rows():
    """
    Test the nestedness_rows function of the Calculations class.
    Check each the value for each pair of rows and for the average Nrow.
    """

    # Test number 1 - 7x4 matrix that appears to be nested
    test_matrix = [[1, 1, 1, 1], 
                   [1, 1, 1, 0],
                   [1, 1, 0, 1],
                   [1, 1, 0, 0],
                   [1, 1, 0, 0],
                   [1, 0, 1, 0],
                   [1, 0, 0, 0]]
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
        val = Calculations(test_matrix, True).nestedness_rows(expected_row_pair)
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
    test_matrix_2 = [[1, 1, 1, 0], 
                     [1, 1, 1, 0],
                     [1, 1, 1, 0],
                     [1, 1, 1, 0],
                     [1, 0, 0, 1],
                     [1, 0, 0, 1],
                     [1, 0, 0, 1]]
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
        val = Calculations(test_matrix_2, True).nestedness_rows(expected_row_pair)
        values.append(val)
        total += val

    assert round(total/21, 1) == expected_val
    assert (values) == [0, 0, 0, 50, 50, 50,
                        0, 0, 50, 50, 50,
                        0, 50, 50, 50,
                        50, 50, 50,
                        0, 0,
                        0]

def test_Calculations_nestedness_cols():
    """
    Test the nestedness_cols function of the Calculations class.
    Check each the value for each pair of columns and for the average Ncol.
    """

    # Test number 1 - 7x4 matrix that appears to be nested
    test_matrix = [[1, 1, 1, 1], 
                   [1, 1, 1, 0],
                   [1, 1, 0, 1],
                   [1, 1, 0, 0],
                   [1, 1, 0, 0],
                   [1, 0, 1, 0],
                   [1, 0, 0, 0]]
    expected_val = 86.1
    total = 0
    values = []
    expected_col_pairs = [(0, 1), (0, 2), (0, 3), 
                          (1, 2), (1, 3), 
                          (2, 3)]
    for expected_col_pair in expected_col_pairs:
        val = Calculations(test_matrix, True).nestedness_cols(expected_col_pair)
        values.append(round(val, 2))
        total += val

    assert round(total/6, 1) == expected_val
    assert (values) == [100, 100, 100,
                        66.67, 100,
                        50]

    # Test number 2 - 7x4 matrix that appears to be modular
    test_matrix_2 = [[1, 1, 1, 0], 
                     [1, 1, 1, 0],
                     [1, 1, 1, 0],
                     [1, 1, 1, 0],
                     [1, 0, 0, 1],
                     [1, 0, 0, 1],
                     [1, 0, 0, 1]]
    expected_val = 50
    total = 0
    values = []
    expected_col_pairs = [(0, 1), (0, 2), (0, 3), 
                          (1, 2), (1, 3), 
                          (2, 3)]
    for expected_col_pair in expected_col_pairs:
        val = Calculations(test_matrix_2, True).nestedness_cols(expected_col_pair)
        values.append(round(val, 2))
        total += val

    assert round(total/6, 1) == expected_val
    assert (values) == [100, 100, 100,
                        0, 0,
                        0]

def test_Calculations_run_parallel():
    """  
    Test the run_parallel function of the Calculations class. 
    Test that nestedness is calculated correctly for the two test matrices.
    """
    # Test number 1 - 7x4 matrix that appears to be nested
    test_matrix = [[1, 1, 1, 1], 
                   [1, 1, 1, 0],
                   [1, 1, 0, 1],
                   [1, 1, 0, 0],
                   [1, 1, 0, 0],
                   [1, 0, 1, 0],
                   [1, 0, 0, 0]]
    expected_val = 82.3
    assert round(Calculations(test_matrix, True).run_parallel(), 1) == expected_val

    # Test number 2 - 7x4 matrix that appears to be modular
    test_matrix_2 = [[1, 1, 1, 0], 
                     [1, 1, 1, 0],
                     [1, 1, 1, 0],
                     [1, 1, 1, 0],
                     [1, 0, 0, 1],
                     [1, 0, 0, 1],
                     [1, 0, 0, 1]]
    expected_val = 39.3
    assert round(Calculations(test_matrix_2, True).run_parallel(), 1) == expected_val



    