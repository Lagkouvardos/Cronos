
How to RUN Chronos
==================



1. Download R Studio (https://www.rstudio.com/) on your computer. 
2. Download the 4 files named: Parameters, Clustering, Transition_Analysis and Modeling.
3. Move the R scripts and the files to one folder all together.
4. Open the R Studio and from File on the upper panel select Open to open the Parameters file. 

At the same time or after you are finished editing your desired parameters on Parameters script, you can open the Modeling file. 

On parameters there are specific comments on what to change and how.

You can select the working directory (folder), which is the directory (folder) where the R scripts and files are saved. The R script has the following format: 
```
#Please set the directory of the script as the working folder (e.g D:/studyname/Data/Chronos/)
#Note: the path is denoted by forward slash "/".
working_directory = "Here_type_the_name_of_the_folder" 
```
##### You should change ONLY what is after the ‘=’ sign. Do not change anything before an ‘=’ sign for the analysis to run. 

Then we change the R working directory (folder) into the one you selected. Then, type the names of the files.
Next to ‘input_otu’ write the name of the tab separated file containing the OTU counts. 
Next to ‘input_meta’ write the name of the file that contains individual sample information and the metadata of the samples.
Next to ‘input_tree’ write the name of the phylogenetic tree constructed from the sequences.

Chronos provides the opportunity to express the clusters that the samples form at a higher taxonomic resolution up to Domain. You can type next to ‘taxonomic_level’ a string denoted by ‘’ signs one of the following: Domain, Phylum, Class, Order or Family. Another option here is to type a number between 1 and 5 respectively. 

Output files and plots produced by Chronos are saved in 2 different directories (folders), the names of which are selected next to dir_with_files and dir_with_plots respectively. 

Finally, the last variable you need to specify is the times the splitting for train and test sets will be performed next to splitting_times. The final accuracy for the train and the test sets are calculated as an average of the splitting times. The optimal percentage for the splitting and overall performance time of the algorithm will also be determined by this number. 

5. Saving the Parameters file by the same name.
6. Open the Modeling file. There all you need to do is select everything (Ctrl+A) and press Run (Ctrl + Enter). 


Chronos scripts description
===========================
Parameters.R: 
______________
This is the file where you select your prefered parameters for Chronos to run. Additionally, files containing the information are specified here and some parameters for the Machine Learning Techniques. 

Clustering.R:
_________________
This script reads all the files for the analysis (! your files will not be altered) and downloads all the required packages in the R environment in order to run the analysis. Samples are assigned to timepoints and the phylogenetic distances between them are calculated via Unifrac metric. Clustering on a dissimilarity matrix calculated for all samples of every timepoint is performed via PAM (Partitioning Around Medoid) method. 

Transition_Analysis.R:
_________________________
Clusters formed in Clustering.R are expressed by their medoid to higher taxonomic resolutions. Samples are assigned to clusters for all timepoints and the transition matrix based on the observations is formed and output as a comma separated file (.csv), which you can open with Excel. In this file, we also check whether the transitions claim the Markovian property.

Modeling.R:
____________________________
This is the script where modeling of the transitions is performed. The analysis is performed in matrices containing all the metadata of the timepoint to predict and the cluster on the timepoint from which prediction is made. Hence, we create matrices, the columns of which are the metadata of the timepoint on which we predict the cluster and the cluster on the timepoint from which prediction is made. Calculations are performed using all the possible combinations of metadata provided, without metadata and using all the metadata available.
We split the dataset into two sets, training set and test set. Training set is used to train/create the model and the test set is where we evaluate it. Training set and test set splits can highly affect the modeling validation metrics, hence we use two independent methods to partition the dataset and repeat the process for a prespecified amount of times (selected on the Parameters.R script). The two methods are Leave One Out (LOO) and Stratified splits. 
LOO is a method where test set is one of the samples and the rest are used to train the model. We apply this partitioning for all the samples on the dataset and output the accuracy of the models as the mean accuracy obtained from the test and train sets.
Stratified splits is a method where a percentage of the samples on the dataset are used to train the model and the rest are used as the test set. The optimal percentage for the splits of all timepoints and all combinations of metadata (or none) are calculated independently, based on the models accuracy. 
This script outputs comma separated files (.csv)  containing all the accuracies achieved and a file with the optimal model as a combination of metadata used. Barplots of the accuracies are output in the plots directory (folder) previously selected on Parameters.R. 
