
![Cronos-logo-v2](https://user-images.githubusercontent.com/8244618/139041662-dc956016-54e9-41b9-bb80-caa91119220a.png)

How to RUN Chronos

Firstly, download R Studio (https://www.rstudio.com/) on your computer. 
Then download the file named Cronos.R. Move the R script and the files to one folder all together.

Open the R Studio and from File on the upper panel select Open to open Cronos.R.
The first lines are dedicated to file and output locations.

You can select the working directory (folder), which is the directory (folder) where the 
R scripts and files are saved. The R script has the following format:

```
# Please set the directory of the script as the working folder (e.g D:/studyname/Data/Chronos/)
# Note: the path is denoted by forward slash "/".
working_directory = "Here_type_the_name_of_the_folder"
```


You should change ONLY what is after the ‘=’ sign. Do not change anything before an ‘=’ sign for the analysis to run. Type the names of the files.
Next to ‘input_otu’ write the name of the tab separated file containing the OTU counts. 
Next to ‘input_meta’ write the name of the file that contains individual sample information and the metadata of the samples.
Next to ‘input_tree’ write the name of the phylogenetic tree constructed from the sequences.

Chronos provides the opportunity to express the clusters that the samples form at a higher taxonomic resolution up to Domain. You can type next to ‘taxonomic_level’ a string denoted by ‘’ signs one of the following: Domain, Phylum, Class, Order or Family. Another option here is to type a number between 1 and 5 respectively. 

Output files and plots produced by Chronos are saved in 4 different directories (folders), within the output directory (folder).

Finally, the last variable you need to specify is the times the splitting for train and test sets will be performed next to splitting_times. The final accuracy for the train and the test sets are calculated as an average of the splitting times. The optimal percentage for the splitting and overall performance time of the algorithm will also be determined by this number. The higher this number is, the more time it will take for Cronos to run.

Now, all you need to do is select everything (Ctrl+A) and press Run (Ctrl + Enter). 


Cronos sections description
Parameters: 
This is the file where you select your prefered parameters for Chronos to run. Additionally, files containing the information are specified here and some parameters for the Machine Learning Techniques. 

Clustering:
This section reads all the files for the analysis (! your files will not be altered) and downloads all the required packages in the R environment in order to run the analysis. Samples are assigned to timepoints and the phylogenetic distances between them are calculated via Unifrac metric. Clustering on a dissimilarity matrix calculated for all samples of every timepoint is performed via PAM (Partitioning Around Medoid) method. 

Transition_Analysis:
Clusters formed in Clustering are expressed by their medoid to higher taxonomic resolutions. Samples are assigned to clusters for all timepoints and the transition matrix based on the observations is formed and output as a comma separated file (.csv), which you can open with Excel. In this file, we also check whether the transitions claim the Markovian property.

Modeling:
This is the section where modeling of the transitions is performed. The analysis is performed in matrices containing all the metadata of the timepoint to predict and the cluster on the timepoint from which prediction is made. Hence, we create matrices, the columns of which are the metadata of the timepoint on which we predict the cluster and the cluster on the timepoint from which prediction is made. Calculations are performed using all the possible combinations of metadata provided, without metadata and using all the metadata available.
We split the dataset into two sets, training set and test set. Training set is used to train/create the model and the test set is where we evaluate it. Training set and test set splits can highly affect the modeling validation metrics, hence we use two independent methods to partition the dataset and repeat the process for a prespecified amount of times (selected on the Parameters section). The two methods are Leave One Out (LOO) and Stratified splits. 
LOO is a method where test set is one of the samples and the rest are used to train the model. We apply this partitioning for all the samples on the dataset and output the accuracy of the models as the mean accuracy obtained from the test and train sets.
Stratified splits is a method where a percentage of the samples on the dataset are used to train the model and the rest are used as the test set. The optimal percentage for the splits of all timepoints and all combinations of metadata (or none) are calculated independently, based on the models accuracy. 
This script outputs comma separated files (.csv)  containing all the accuracies achieved and a file with the optimal model as a combination of metadata used. Barplots of the accuracies are output in the Accurcies directory (folder).  
