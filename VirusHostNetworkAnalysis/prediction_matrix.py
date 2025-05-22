import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

class PredictionMatrix:
    """ Class to make a matrix from the input file. Can make square or rectangle matrices. Can also be used to make probabilty matrix.
    
    Args:
    vhip_predictions_file (str): Path to the input file containing virus-host interaction predictions.
    
    """

    def __init__(self, vhip_predictions_file:str, probability:bool = False):
        #check if this file path exists
        if not os.path.exists(vhip_predictions_file):
            raise FileNotFoundError(f"The file {vhip_predictions_file} does not exist.")
        else:
            self.virus_host = pd.read_csv(vhip_predictions_file, sep='\t')
            self.error_check()
            self.title = vhip_predictions_file.split('/')[-1].split('.')[0]
            # If probability is True, then the matrix will be a probability matrix
            self.probability = probability

    def error_check(self):
        """ Check if the input file is in the correct format. """
        # Check if the file has the three columns needed for visualization
        # The files need pairs, Predicitons, and either InfProbabilities OR Scores
        required_columns = ['pairs', 'Predictions']
        if not all(col in self.virus_host.columns for col in required_columns):
            raise ValueError(f"The input file must contain the following columns: {', '.join(required_columns)}")
        # Check if the file has either InfProbabilities or Scores
        if not any(col in self.virus_host.columns for col in ['InfProbabilities', 'Scores']):
            raise ValueError("The input file must contain either the InfProbabilities or Scores column.")
        # Check if the pairs column is in the correct format
        if not all(self.virus_host['pairs'].str.contains(':')):
            raise ValueError("The pairs column is not in the correct format. It should be in the format 'virus:host'")
    
    def plot_scores(self):
        """ Plot the scores of the predictions. """
        plt.figure(figsize=(10, 6))
        # column name is either InfProbabilities or Scores
        if 'InfProbabilities' in self.virus_host.columns:
            sns.histplot(self.virus_host['InfProbabilities'], kde=False)
        else:
            sns.histplot(self.virus_host['Scores'], kde=False)
        plt.title(f"Distribution of Scores for {self.title}")
        plt.xlabel("Scores")
        plt.ylabel("Frequency")
        plt.show()

    def get_unique_virus_host(self):
        """Find the unique virus names and host names in the dataset. The first column contains a virus name and host name, which are separated by a colon."""
        # The unique viruses will make of the rows.
        self.rows = self.virus_host['pairs'].str.split(':').str[0].unique()
        # The unique hosts will make up the columns.
        self.columns = self.virus_host['pairs'].str.split(':').str[1].unique()

    def initialize_matrix(self):
        """Create a matrix full of zeros and make a list of rows and column labels. """
        self.virus_host_array = np.zeros((len(self.rows), len(self.columns)), dtype= bool if self.probability is False else float)
    
    def fill_matrix(self):
        """ Fill the prediction matrix with 1s and 0s or fill the probability matrix with the InfProbabilities from the dataset."""
        if self.probability is False:
            # Only loops through the subset of data where predicitons == 1
            for _, row in self.virus_host[self.virus_host['Predictions'] == 1].iterrows():
                virus, host = row['pairs'].split(':')
                # Find the index of the virus (row) and host (col) in the matrix
                virus_index = np.where(self.rows == virus)[0][0]
                host_index = np.where(self.columns == host)[0][0]
                # Fill the matrix with 1 if matrix_type == 'predictions'
                self.virus_host_array[virus_index][host_index] = 1
        else:
            for _, row in self.virus_host.iterrows():
                virus, host = row['pairs'].split(':')
                # Find the index of the virus (row) and host (col) in the matrix
                virus_index = np.where(self.rows == virus)[0][0]
                host_index = np.where(self.columns == host)[0][0]
                # Fill the matrix with InfProbabilities or Scores depending on the column name
                if 'InfProbabilities' in row:
                    self.virus_host_array[virus_index][host_index] = row['InfProbabilities']
                elif 'Scores' in row:
                    self.virus_host_array[virus_index][host_index] = row['Scores']

    def make_rectangular_matrix(self):
        """ Call all functions to make a rectangular matrix
        
        Returns:
            np.ndarray: The filled matrix with the virus-host interactions.
        """
        self.get_unique_virus_host()
        self.initialize_matrix()
        self.fill_matrix()
        return self.virus_host_array

    def save_matrix(self):
        """ Save the matrix to a csv file. """
        np.savetxt('virus_host.csv', self.virus_host_array, delimiter=',')
        print('virus_host.csv saved')