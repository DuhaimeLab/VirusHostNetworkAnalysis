The data will look like:
* Results of VHIP 2.0 in .tsv form
    * Column 1 will have pairs data, where the side of the colon is the virus and the right side of the colon is the host.
    * Another column for Predictions will be used to connect the nodes for virus and host
    
 * Link to the data: https://github.com/ecampau/VirusHostNetworkAnalysis/tree/main/Sample_Input
 
 * Validation stats:
 	* The data is the result output from VHIP 2.0, so we already know the columns will be in the same order each time. 
	* Random models, such as ER and Configuration models, will be constructed and used as a null model to validate the results of the virus-host network model. See null_models.py for these models.
	* Each set of data has more than one unique host and one unique virus