
![Cronos-logo-v2](https://user-images.githubusercontent.com/8244618/139041662-dc956016-54e9-41b9-bb80-caa91119220a.png)

## Introduction
Microbial time-series analysis, typically, examines the abundances of individual taxa over time and attempts to assign etiology to observed patterns. Here, we suggest an alternative approach to the analysis of microbial time-series. We attempt to track the transition of communities, rather than individual taxa. Cronos detects the intrinsic microbial profile clusters on all time points, describes them in terms of composition, and records the transitions between them. Cluster assignments, combined with the provided metadata, are used to model the transitions and predict samples’ fate under various effects. 

## Requirements and Installation

### Input Files
The user has to provide three mandatory files

* An OTUs or ASVs abundance table.
* A phylogenetic tree that corresponds to the OTUs or ASVs of the abundance table.
* A mapping file that contains the labels of the samples and the information about the Time points.

The pre-processing of the raw sequencing data through the [IMNGS platform](www.imngs.org/ "IMNGS site") and [Rhea pipeline](https://github.com/Lagkouvardos/Rhea/ "Rhea repository") can provide the user with the necessary input files for the Cronos.

### How to RUN Cronos

In order to execute the script, it is required the [R language](https://www.r-project.org/ "R download site") to be locally installed. The use of the [R-Studio](https://www.rstudio.com/products/rstudio-desktop/ "R-studio download site") will simplify the procedure even for non-experienced users. During the first execution of the program, all the necessary packages will be installed, so at least for this first run, a stable internet connection is required.

Cronos does not require any installation process.

In order to run Cronos the user has to follow the next steps:
*	First, download the Cronos project from the Github repository.
*	Then decompress the files into the desired destination.
* Open the Cronos.R script
* Set the working directory where the files are saved
* Fill out the required parameters as they detailed described in the script
* Run the script

## Script structure

Each script contains: 
* The Script Overview section
* The Initialization section
* The Main section

The general rule that applies is that the lines concerning the Overview and Initialization sections start with a hash sign followed by a back-tick (#`). These lines aim to inform the user about the details of the program and the actions that should be followed.

You should change ONLY what is after the ‘=’ sign in the Initialization section. Do not change anything before an ‘=’ sign. 


## Outputs

All the Output tables and plots produced by Cronos are saved in a folder in the same dicetrory with the original input files.




 
## Citation

Litos, A., Intze, E., Pavlidis, P., & Lagkouvardos, I. (2022). Cronos: A Machine Learning Pipeline for Description and Predictive Modeling of Microbial Communities Over Time. Frontiers in Bioinformatics. https://doi.org/10.3389/fbinf.2022.866902
