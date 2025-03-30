The data will look like:
* Results of VHIP 2.0 in .tsv form
    * Column 1 will have pairs data, where the side of the colon is the virus and the right side of the colon is the host.
    * Another column for Predictions will be used to connect the nodes for virus and host.
	* A third column, Probability, will be used to make heatmaps.
    
 * Link to the data: https://github.com/ecampau/VirusHostNetworkAnalysis/tree/main/Sample_Input
 * The sample dataset Aug4_predictions.tsv can be used to follow the tutorial and is a subset of the total provided input data. After running with the sample subset, the package will be run with the other files in the Sample_Input folder, as well as other datasets not yet uploaded.
 
 * Validation stats:
 	* The data is the result output from VHIP 2.0, so we already know the columns will be in the same order each time. To double check this, there are some functions to check that the columns that are needed are included and in the correct format. 
	* Random models, such as ER and Configuration models, will be constructed and used as a null model to validate the results of the virus-host network model. See null_models.py for these models.
	* Each set of data has more than one unique host and one unique virus.