##########################################################################################
######## IN THIS SECTION THE USER CAN SET PARAMETERS AND FILE LOCATIONS ##################
##########################################################################################

########### PLEASE FOLLOW THE INSTRUCTIONS CAREFULLY #####################################

#' Please set the directory of the script as the working folder (e.g D:/studyname/Data/Chronos/)
#' Note: the path is denoted by forward slash "/"
setwd("~/Working_Chronos/")  #<--- CHANGE ACCORDINGLY !!!

#' Please give the file name of the normalized OTU-table without taxonomic classification
input_otu = "OTUs_Table.csv"              #<--- CHANGE ACCORDINGLY !!!
#' Please give the name of the meta-file that contains individual sample information
input_meta = "Map_OTUs"                #<--- CHANGE ACCORDINGLY !!!
#' Please give the name of the phylogenetic tree constructed from the OTU sequences
input_tree = "OTUs-NJTree.tre"                   #<--- CHANGE ACCORDINGLY !!!

# Please select the taxon in which the samples will be analyzed
# Either type it in ' ' e.g. 'Order' or select a number between 1 and 5, where:
# 1: Domain,
# 2: Phylum
# 3: Class
# 4: Order
# 5: Family

taxonomic_level=3  #### <---- CHANGE ACCORDINGLY


# Please select method with which the optimal number of clusters will be selected
# Could be either Drop or Highest
method='Drop'

######################### END OF SECTION #################################################


#################### NO CHANGE IS NEEDED BEYOND HERE #####################################

#### You can select everything (Ctrl+A) and press run (Ctrl+Enter)  


##########################################################################################
############################ READING THE FILES ###########################################
##########################################################################################

meta_file <- read.table (file = input_meta, check.names = FALSE, header = TRUE, dec = ".", sep = ",", row.names = 1, comment.char = "")
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
  return (taxa)
}
taxa_matrix<- taxa_matrix(otus_taxonomic,otu_file)

######################### END OF SECTION #################################################

##########################################################################################
################# CONVERT FILES TO DESIRABLE FORMAT ######################################
##########################################################################################


# keep only those rows that appear in the mapping file
otu_file <- otu_file[,rownames(meta_file)]

# OTU-table and mapping file should have the same order and number of sample names
# Order the OTU-table by sample names (ascending)
otu_file <- otu_file[,order(names(otu_file))]
# Transpose OTU-table and convert format to a data frame
otu_file <- data.frame(t(otu_file))


######################### END OF SECTION #################################################

##########################################################################################
############### CHECKING FOR AND INSTALLING PACKAGES  ####################################
##########################################################################################

packages <-c("ade4","GUniFrac","phangorn","cluster","fpc") 
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

### EDO NA MPEI H SYNARTISI POY UA TO SPAEI, OTAN BRV NA APOTHIKEYV DIAFORETIKOY MEGETHOUS
### PRAGMATA KAPOU


######################### END OF SECTION #################################################

##########################################################################################
############# CALCULATE ALL UNIFRAC DISTANCES OF SAMPLES PER TIMEPOINT ###################
##########################################################################################
# Calculate the UniFrac distance matrix for comparing microbial communities
unifracs <- GUniFrac(otu.tab = otu_file ,tree = rooted_tree, alpha = c(0.0,0.5,1.0))$unifracs

# Weight on abundant lineages so the distance is not dominated by highly abundant lineages with 0.5 having the best power
unifract_dist <- unifracs[, , "d_0.5"]
unifract_dist
######################### END OF SECTION #################################################


##########################################################################################
############# CLUSTERING #################################################################
##########################################################################################

k=4

################ PRACTICE ############################

kentra<- sample(x = row.names(unifract_dist), size = k)

clustering_function <-function(unifract_dist,k,kentra){
  
  
  distances <- function(unifract_dist,kentra,k) {
    
    apostaseis<- matrix(0,nrow = nrow(unifract_dist), ncol = k)
    for (i in 1:nrow(unifract_dist)){
      for (j in 1:k){
        apostaseis[i,j]= unifract_dist[i,kentra[j]]
      }
      apostaseis[i,]=apostaseis[i,]==min(apostaseis[i,])
    }
    return (apostaseis)
  }
  
  medoids<- function(unifract_dist,apostaseis,kentra){
    for (i in 1:k){
      kentra[i] <- row.names(unifract_dist[apostaseis[,i]==1,apostaseis[,i]==1])[which.min(apply(X = unifract_dist[apostaseis[,i]==1,apostaseis[,i]==1], FUN = mean,MARGIN=2))]
      }
    return (kentra)
  }
  
  matequal <- function(x, y){
    return (is.matrix(x) && is.matrix(y) && dim(x) == dim(y) && all(x == y))
  }
  
  apostaseis <- matrix(1,nrow = nrow(unifract_dist), ncol = k)
  palies_apostaseis <- matrix(0,nrow = nrow(unifract_dist), ncol = k)
  palia_kentra<- c()
  while (!(matequal(apostaseis,palies_apostaseis) && all(kentra==palia_kentra))) {
    palia_kentra<- kentra
    palies_apostaseis<- apostaseis
    apostaseis<- distances(unifract_dist = unifract_dist, kentra = kentra,k = k)
    kentra<- medoids(unifract_dist = unifract_dist, apostaseis = apostaseis, kentra = kentra)
  }
  return (apply(apostaseis, FUN = which.max, MARGIN = 1))
}



optimal_k<- function(unifract_dist, method){
  calinski_harabasz_values <-c()
  for (k in 3:12){
    kentra<- sample(x = row.names(unifract_dist), size = k)
    calinski_harabasz_values= c(calinski_harabasz_values,(calinhara(x = unifract_dist,cn = k, clustering = clustering_function(unifract_dist = unifract_dist,k = k,kentra = kentra ))))
  }
  print (diff(calinski_harabasz_values))
  if (method=='Drop'){
    kalitero<- which.min(diff(calinski_harabasz_values))+1
  }
  else if (method=='Highest'){
    kalitero<- which.max(calinski_harabasz_values)
  }
  return (kalitero)
}

k=optimal_k(unifract_dist = unifract_dist, method = method)
k
clusters <- clustering_function(unifract_dist = unifract_dist,k = k, kentra = sample(x = row.names(unifract_dist), size = k))


samples_on_clusters<-list()
for (i in 1:k){
  samples_on_clusters[[i]] <- as.vector(row.names(otu_file)[clusters==k])
}
samples_on_clusters[2]
