##########################################################################################
############ IN THIS FILE YOU CAN SET PARAMETERS AND FILE LOCATIONS ######################
##########################################################################################

########### PLEASE FOLLOW THE INSTRUCTIONS CAREFULLY #####################################

#' Please set the directory of the script as the working folder (e.g D:/studyname/Data/Chronos/)
#' Note: the path is denoted by forward slash "/"
setwd("~/Working_Chronos/Chronos_almost")              #<--- CHANGE ACCORDINGLY !!!

#' Please give the file name of the normalized OTU-table without taxonomic classification
input_otu = "SOTUs-Table.tab"           #<--- CHANGE ACCORDINGLY !!!
#' Please give the name of the meta-file that contains individual sample information
input_meta = "Mapping_File_Inf_St.csv"                #<--- CHANGE ACCORDINGLY !!!
#' Please give the name of the phylogenetic tree constructed from the OTU sequences
input_tree = "SOTUs-NJTree.tre"         #<--- CHANGE ACCORDINGLY !!!

# Please specify if the file contains both adult and infant data: 
# If it contains adult data specify the name of the column (i.e. timepoint name) where they are saved
adult_timepoint_name = 'MM'              #<--- CHANGE ACCORDINGLY


# Please select the taxon in which the samples will be analyzed
# Either type it in ' ' e.g. 'Order' or select a number between 1 and 5, where:
# 1: Domain,
# 2: Phylum
# 3: Class
# 4: Order
# 5: Family

taxonomic_level='Family'                          # <---- CHANGE ACCORDINGLY


# Please select method for cluster representation. It can be either mean, or median.
# Cluster representation via its mean calculates the mean abundancy of every taxon in the 
# cluster and return the mean. It is preferable in larger datasets and creates a composition
# not found in any of the samples on the cluster.
# Cluster representation via its median calculates the median of the samples belonging to
# the cluster. It is preferable in smaller datasets and returns a illustative composition
# that belongs to a sample.
representation_method = 'median'

# Please select clustering method. It could be one of hierarchical or 
# PAM (Partition Around Medoids). If your dataset is small, we recommend hierarchical
# else, we recommend PAM. PAM is also recommended for expected well-separated datasets
clustering_method = 'Hierarchical'            # <---- CHANGE ACCORDINGLY


# Please write the names of the 2 new directories, in which the output should be saved at:
# Firstly, the one for the transition plots
dir_with_plots= 'Plots'     # <---- CHANGE ACCORDINGLY
# Secondly, the one with the files
dir_with_files= 'Chronos_output_files' # <---- CHANGE ACCORDINGLY