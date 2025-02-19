from VirusHostNetworkAnalysis import prediction_matrix

def test_square_matrix():
    # Test the unique identifiers
    test_matrix = prediction_matrix.PredictionMatrix('Sample_Input/test_predictions.tsv')
    test_matrix.get_unique_virus_host()
    assert len(test_matrix.unique_viruses) == 2

def test_filling_correctly():
    # Test the matrix is filled correctly
    test_matrix = prediction_matrix.PredictionMatrix('Sample_Input/test_predictions.tsv')
    # fill matrix for predictions and test that it only has 1s and 0s
    test_matrix.get_unique_virus_host()
    test_matrix.initialize_matrix('prediction')
    test_matrix.fill_matrix('prediction')
    # assert that matrix contains only 1s and 0s
    assert test_matrix.virus_host_array.all() == 1 or test_matrix.virus_host_array.all() == 0

def test_sorting():
    # Test the matrix is sorted correctly
    test_matrix = prediction_matrix.PredictionMatrix('Sample_Input/test_predictions.tsv')
    test_matrix.get_unique_virus_host()
    test_matrix.initialize_matrix('prediction')
    test_matrix.fill_matrix('prediction')
    test_matrix.sort_matrix()
    # assert that the matrix is sorted correctly
    assert test_matrix.virus_host_array[0][0] == 1
    assert test_matrix.virus_host_array[-1][-1] == 0
