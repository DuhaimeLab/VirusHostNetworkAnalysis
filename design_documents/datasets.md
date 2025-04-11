The data will look like:
* Results of VHIP 2.0 in .tsv form
    * Column 1 will have pairs data, where the side of the colon is the virus and the right side of the colon is the host.
    * Another column for Predictions will be used to connect the nodes for virus and host.
	* A third column, Probability, will be used to make heatmaps.
    
 * Link to the data: https://github.com/ecampau/VirusHostNetworkAnalysis/tree/main/Sample_Input
 * Example dataset for using the tool:
	* A small set of data was created that can be used to learn the tool and more easily visualize the ouput. This dataset is called data_subset.tsv and has just 6 viruses and 3 hosts. Viruses and hosts were taken from Sep29_predictions.tsv randomly. This was mostly done to create graphs for the ER and Configuration models that have few nodes and edges and allow the user to see changes. This was run in the tutorial, followed by a larger subset of data.
 	* The sample dataset Aug4_predictions.tsv can also be used to follow the tutorial and is a subset of the total provided input data (larger subset than data_subset.tsv).  In the tutorial for the true data, this package will be run with the other files in the Sample_Input folder, and in the future it will include other datasets not yet uploaded. The file contains 108,937 rows of virus-host pairs. After running pre-processing steps and organiing into a mtrix of unique values, the matrix size if 6408 by 17. This data was chosen as an example dataset because the resulting matrix shows a clear and interesting pattern that might not be present or as visible in the other data subsets. Additionally, this data has been cleaned and we know there are no errors, missing values, or duplicates included.


 * Real dataset for answering a biological question using the tool.
	* The real dataset is the full data included in the Sample_Input folder. This data is a result of VHIP being run on samples from Lake Erie. The samples were originally collected to study cyanobacterial harmful algae blooms in the water. Harmful algal blooms are the rapid growth of algae or cyanobacteria in water that can harm people, animals, or the environment. Cyanobacterial harmful algae blooms (cHABs) are specific to freshwater and are characterized by their blue-green appearance. The data from Lake Erie was collected by AJ Wing and the analysis was run by Eric Bastien and Evelyn Faust. VHIP 2.0 aims to predict which viruses interact which hosts. The output is this prediction data from the site at varying time points. Moving forward, we hope to get access to data for a longer span of time and run this analysis on each of those additional time points. This will allow us to address another question: "How do patterns in predicted virus-host interactions change over time?". The expected results can be found in VHIP_real_data_tutorial.ipynb.

 * Validation stats:
 	* The data is the result output from VHIP 2.0, so we already know the columns will be in the same order each time. To double check this, there are some functions to check that the columns that are needed are included and in the correct format. 
	* Random models, such as ER and Configuration models, will be constructed and used as a null model to validate the results of the virus-host network model. See null_models.py for these models.
	* Each set of data has more than one unique host and one unique virus.
