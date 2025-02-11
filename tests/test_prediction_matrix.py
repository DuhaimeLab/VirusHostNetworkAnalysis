import os
import sys
from VirusHostNetworkAnalysis import prediction_matrix

def test_square_matrix():
    # Test the unique identifiers
    test_matrix = prediction_matrix.PredictionMatrix('Sample_Input/test_predictions.tsv')
    test_matrix.get_unique_virus_host()
    assert len(test_matrix.unique_viruses) == 2