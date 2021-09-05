##########################################################################################
######## IN THIS SECTION THE USER CAN SET PARAMETERS AND FILE LOCATIONS ##################
##########################################################################################

########### PLEASE FOLLOW THE INSTRUCTIONS CAREFULLY #####################################

#' Please set the directory of the script as the working folder (e.g D:/studyname/Data/Chronos/)
#' Note: the path is denoted by forward slash "/"
setwd("~/Working_Chronos/Chronos_almost")              #<--- CHANGE ACCORDINGLY !!!

#' Please give the file name of the normalized OTU-table without taxonomic classification
input_otu = "OTUs-TableIS.tab"           #<--- CHANGE ACCORDINGLY !!!
#' Please give the name of the meta-file that contains individual sample information
input_meta = "Meta_Inf_Stud.tab"                #<--- CHANGE ACCORDINGLY !!!
#' Please give the name of the phylogenetic tree constructed from the OTU sequences
input_tree = "OTUs-NJTree.treIS"         #<--- CHANGE ACCORDINGLY !!!

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

taxonomic_level=4                          # <---- CHANGE ACCORDINGLY


# Please select method with which the optimal number of clusters will be selected
# Could be either Drop or Highest
clustering_method='Highest'                # <---- CHANGE ACCORDINGLY

# Please select method to declare the transition matrix
# It can be: mle (Maximum Likelihood Estimation), map (Maximum a posteriori) or  bootstrap
markov_method= 'mle'

# Please write the names of the 2 new directories, in which the output should be saved at:
# Firstly, the one for the transition plots
dir_with_plots= 'Transition_Plots'     # <---- CHANGE ACCORDINGLY
# Secondly, the one with the files
dir_with_files= 'Chronos_output_files' # <---- CHANGE ACCORDINGLY


######################### END OF SECTION #################################################


#################### NO CHANGE IS NEEDED BEYOND HERE #####################################

#### You can select everything (Ctrl+A) and press run (Ctrl+Enter)  


##########################################################################################
############################ READING THE FILES ###########################################
##########################################################################################

meta_file <- read.table (file = input_meta, check.names = FALSE, header = TRUE, dec = ".", sep = "\t", row.names = 1, comment.char = "", stringsAsFactors = F)
# Clean table from empty lines
meta_file <- data.frame(meta_file[!apply(is.na(meta_file) | meta_file=="",1,all),])
# Order the mapping file by sample names (ascending)
meta_file <- data.frame(meta_file[order(row.names(meta_file)),])
# Load the tab-delimited file containing the values to be analyzed (samples names in the first column)
otu_file <- read.csv(file = input_otu,sep = '\t',row.names = 1,header = T, stringsAsFactors = F)
# Clean table from empty lines
otu_file <- otu_file[!apply(is.na(otu_file) | otu_file =="",1,all),]

######################### END OF SECTION #################################################

##########################################################################################
### THIS SECTION EXPRESSES THE TAXONOMIC PROFILE OF EACH SAMPLE AT THE SELECTED LEVEL ####
##########################################################################################

taxonomy_of_sample<- function(taxonomic_level,otu_file){
  
  taxonomic_levels<- c('Domain','Phylum','Class','Order','Family','Genus')
  taxonomies<-c()
  
  if (typeof(taxonomic_level)=='double'){
    taxonomic_level= taxonomic_levels[taxonomic_level]
}
  for (i in row.names(otu_file)){
    otu_file[i,'taxonomy']<-unlist(strsplit(otu_file[i,'taxonomy'], split = ';'))[which(taxonomic_levels==taxonomic_level)]

  }
  otu_file[,'taxonomy']
  return (otu_file)
}

otus_taxonomic = taxonomy_of_sample(taxonomic_level,otu_file)


taxa_matrix<-function(otus_taxonomic,otu_file){
  taxa_selected<- unique(otus_taxonomic[,'taxonomy'])
  
  taxa<- matrix(0,nrow = length(taxa_selected),ncol = ncol(otus_taxonomic)-1)
  colnames(taxa)<- head(colnames(otu_file),-1)
  row.names(taxa)<- taxa_selected
  
  for (i in 1:length(taxa_selected)){
    taxa[i,] = colSums(otus_taxonomic[otus_taxonomic['taxonomy']==taxa_selected[i],1:ncol(otus_taxonomic)-1])
  }
  return (prop.table(x = taxa,margin = 2))
}
taxa_matrix<- taxa_matrix(otus_taxonomic,otu_file)
taxa_matrix <- data.frame(t(taxa_matrix))

######################### END OF SECTION #################################################

##########################################################################################
################# CONVERT FILES TO DESIRABLE FORMAT ######################################
##########################################################################################
kati <- c()
kati_allo <-c()
for (i in rownames(meta_file)){
  kati<- append(kati, paste(ifelse(nchar(i)==5,yes = 'X',no = 'X0'), i, sep = ""))
  kati_allo <- append(kati_allo, paste(ifelse(nchar(meta_file[i,'Sample'])==2,yes = "X0",no = "X"),meta_file[i,'Sample'], sep = "")) 
}
rownames(meta_file) <- kati
meta_file[,'Sample'] <- kati_allo
meta_file
# keep only those rows that appear in the mapping file
otu_file <- otu_file[,rownames(meta_file)]
# OTU-table and mapping file should have the same order and number of sample names
# Order the OTU-table by sample names (ascending)
otu_file <- otu_file[,order(names(otu_file))]
# Transpose OTU-table and convert format to a data frame
otu_file <- data.frame(t(otu_file))
# keep only those that appear in the otu file
meta_file <- meta_file[rownames(otu_file),]

######################### END OF SECTION #################################################

##########################################################################################
############### CHECKING FOR AND INSTALLING PACKAGES  ####################################
##########################################################################################

packages <-c("ade4","GUniFrac","phangorn","cluster","fpc","markovchain") 
# Function to check whether the package is installed
InsPack <- function(pack)
{
  if ((pack %in% installed.packages()) == FALSE) {
    install.packages(pack,repos ="http://cloud.r-project.org/")
  } 
}

# Applying the installation on the list of packages
lapply(packages, InsPack)
# Make the libraries
lib <- lapply(packages, require, character.only = TRUE)
# Check if it was possible to install all required libraries
flag <- all(as.logical(lib))


# Load the phylogenetic tree calculated from the OTU sequences 
tree_file <- read.tree(input_tree)
# Root the OTU tree at midpoint 
rooted_tree <- midpoint(tree_file)

######################### END OF SECTION #################################################

##########################################################################################
############# DIVIDE THE DATASET INTO DIFFERENT TIMEPOINTS ###############################
##########################################################################################

timepoint_collection <- function(otu_file,meta_file){
  timepoints = unique(meta_file[,'Timepoint'])
  timepoint_list<- list()
  for (i in timepoints){
    timepoint_list[[as.character(i)]] <- otu_file[meta_file[,'Timepoint']==i,]
  }
  timepoint_list
  return (timepoint_list)
}
timepoint_list <- timepoint_collection(otu_file,meta_file = meta_file)


# Create a matrix to save the sample clusters by timepoint
samples_on_clusters<-matrix(0, ncol = length(timepoint_list), nrow= nrow(unique(meta_file['Sample'])))
row.names(samples_on_clusters)<- unique(meta_file[,'Sample'])
colnames(samples_on_clusters)<- names(timepoint_list)

######################### END OF SECTION #################################################

##########################################################################################
############# CALCULATE ALL UNIFRAC DISTANCES OF SAMPLES PER TIMEPOINT ###################
##########################################################################################
# Calculate the UniFrac distance matrix for comparing microbial communities
for (name in names(timepoint_list)){

unifracs <- GUniFrac(otu.tab = timepoint_list[[name]] ,tree = rooted_tree, alpha = c(0.0,0.5,1.0))$unifracs

# Weight on abundant lineages so the distance is not dominated by highly abundant lineages with 0.5 having the best power
unifract_dist <- unifracs[, , "d_0.5"]

######################### END OF SECTION #################################################

##########################################################################################
############# CLUSTERING #################################################################
##########################################################################################

k=2
kentra<- sample(x = row.names(unifract_dist), size = k)


clustering_function <-function(unifract_dist,k,kentra){
  
  lis <- list()
  distances <- function(unifract_dist,kentra,k) {
    lista<- list()
    egine_sosta<- T
    apostaseis<- matrix(0,nrow = nrow(unifract_dist), ncol = k)
    for (i in 1:nrow(unifract_dist)){
      for (j in 1:k){
        apostaseis[i,j]= unifract_dist[i,kentra[j]]
      }
      apostaseis[i,]=apostaseis[i,]==min(apostaseis[i,])
    }
    if (any(colSums(apostaseis)==1)){
      egine_sosta<- F
      kentra<- sample(x = row.names(unifract_dist), size = k)
      return (distances(unifract_dist,kentra ,k))
    }
    
    lista[[1]]<- apostaseis
    lista[[2]]<- egine_sosta
    return (lista)
  }
  
  medoids<- function(unifract_dist,apostaseis,kentra,k){
    for (i in 1:k){
      kentra[i] <- row.names(unifract_dist[apostaseis[,i]==1,apostaseis[,i]==1])[which.min(apply(X = unifract_dist[apostaseis[,i]==1,apostaseis[,i]==1], FUN = mean,MARGIN=1))]
      }
    return (kentra)
  }
  
  matequal <- function(x, y){
    return (dim(x) == dim(y) && all(x == y))
  }
  
  apostaseis <- matrix(1,nrow = nrow(unifract_dist), ncol = k)
  palies_apostaseis <- matrix(0,nrow = nrow(unifract_dist), ncol = k)
  palia_kentra<- c()
  while (!(matequal(apostaseis,palies_apostaseis) && all(kentra==palia_kentra))) {
    palia_kentra<- kentra
    palies_apostaseis<- apostaseis
    apotelesma<- distances(unifract_dist = unifract_dist, kentra = kentra,k = k)
    apostaseis <- apotelesma[[1]]
    egine_sosta <- apotelesma[[2]]
    kentra<- medoids(unifract_dist = unifract_dist, apostaseis = apostaseis, kentra = kentra, k = k)
  }
  
  lis[[1]]<- apply(apostaseis, FUN = which.max, MARGIN = 1)
  lis[[2]]<- egine_sosta
  
  return (lis)
}



optimal_k<- function(unifract_dist, method){
  calinski_harabasz_values <-c()
  for (k in 2:8){
    kentra<- sample(x = row.names(unifract_dist), size = k)
    apotelesma <- clustering_function(unifract_dist = unifract_dist,k = k,kentra = kentra )
    clustering <-apotelesma[[1]]
    egine_sosta<- apotelesma[[2]]
    calinski_harabasz_values= c(calinski_harabasz_values,(calinhara(x = unifract_dist,cn = k, clustering = clustering)))
    if (egine_sosta == F){
      break
    }
  }
  if (method=='Drop'){
    kalitero<- which.min(diff(calinski_harabasz_values))+1
  }
  else if (method=='Highest'){
    kalitero<- which.max(calinski_harabasz_values)
  }
  return (kalitero)
}

k=optimal_k(unifract_dist = unifract_dist, method = clustering_method)
clusters <- clustering_function(unifract_dist = unifract_dist,k = k, kentra = sample(x = row.names(unifract_dist), size = k))[[1]]

samples_on_clusters[meta_file[row.names(unifract_dist),'Sample'],name]<- clusters
}

######################### END OF SECTION #################################################

##########################################################################################
######################## DECLARE TRANSITION MATRIX  ######################################
##########################################################################################

# Exclude the infants of the dataset
infants<- samples_on_clusters[samples_on_clusters[,adult_timepoint_name]==0,1:ncol(samples_on_clusters)-1]

# Create the markovian transition matrix for each combination of timepoints and save them

transition_matrices <- list()
for (i in colnames(infants)){
  for (j in colnames(infants)){
    if (as.numeric(i) < as.numeric(j)){
      markovestimation <- markovchainFit(as.character(infants[,c(i,j)]), method = markov_method, byrow = T)
      transition_matrices[[paste(i,j,sep = '_')]] <- t(markovestimation$estimate)
    }
  }
}
######################### END OF SECTION #################################################

##########################################################################################
######################## SPECIFY TAXA ON CLUSTERS   ######################################
##########################################################################################


taxa_per_cluster <- function(taxa_matrix,samples_on_clusters,timepoint_list){
  taxa_clusters <- list()
  for (i in names(timepoint_list)){
    for (j in 1:max(samples_on_clusters[,i])){
      taxa_clusters[[paste(as.character(i),as.character(j),sep = ' cluster ')]] <- as.matrix(x = apply(X = taxa_matrix[paste(rownames(samples_on_clusters[samples_on_clusters[,i]==j,]),i,sep = ''),], MARGIN = 2, FUN = mean),decreasing = T)      
      }
  }
  clustering_taxa <- matrix(0,  nrow = ncol(taxa_matrix), ncol = length(taxa_clusters))
  row.names(clustering_taxa)<- colnames(taxa_matrix)
  colnames(clustering_taxa)<- names(taxa_clusters)
  for (i in colnames(clustering_taxa)){
    for (j in rownames(clustering_taxa)){
      clustering_taxa[j,i] <- taxa_clusters[[i]][j,]
    }
  }
  return (clustering_taxa)
}


taxa_clusters <- taxa_per_cluster(taxa_matrix = taxa_matrix, samples_on_clusters = samples_on_clusters, timepoint_list = timepoint_list)

######################### END OF SECTION #################################################

##########################################################################################
################ PLOT TRANSITION PROBABILITIES ON DIRECTORY  #############################
##########################################################################################

dir.create(dir_with_plots)

for (name in names(transition_matrices)){
  jpeg(filename =paste('Transition_Plots', paste(paste(unlist(strsplit(name,split = '_'))[1],unlist(strsplit(name,split = '_'))[2], sep = ' to '), 'timepoints', sep = ' '),sep = '/'))
  plot(transition_matrices[[name]], name = name)
  dev.off()
  
}

######################### END OF SECTION #################################################


##########################################################################################
################ WRITE TAB DELIMITED FILES WITH THE OUTPUTS  #############################
##########################################################################################


dir.create(dir_with_files)


colnames(taxa_clusters)<-lapply(X = colnames(taxa_clusters), FUN = function(x){paste('Timepoint',x,sep = ' ')})
row.names(samples_on_clusters) <- lapply(X = rownames(samples_on_clusters), FUN = function(x){substr(x,start=2,stop= 4)})

write.csv(x = samples_on_clusters, file = paste(dir_with_files, "Samples_in_Timepoint-specific_Clusters.csv",sep = '/'), row.names = T)
write.csv(x = taxa_clusters,       file = paste(dir_with_files, "Taxonomic_profile_of_clusters.csv", sep = '/'), row.names = T)

######################### END OF SECTION #################################################

############## PRACTICE ############################

############## COMMENTS ############################
# Καταλήγω με πίνακα που θα έχει τα διαφορετικά clusters ανά δείγμα και ως features ????
# 27 δείγματα είναι πολύ λίγα για random forest, ειδικά χωρίς features.
