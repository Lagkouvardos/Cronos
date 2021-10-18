#################### NO CHANGE IS NEEDED HERE ############################################

#### You can select everything (Ctrl+A) and press run (Ctrl+Enter)  

source('Parameters.R')
dir.create(dir_with_plots, showWarnings = F)


############################ Reading the files ###########################################
meta_file <- read.table (file = input_meta, check.names = FALSE, header = TRUE, dec = ".", sep = ",", row.names = 1, comment.char = "", stringsAsFactors = F)
# Clean table from empty lines
meta_file <- data.frame(meta_file[!apply(is.na(meta_file) | meta_file=="",1,all),])
# Order the mapping file by sample names (ascending)
meta_file <- data.frame(meta_file[order(row.names(meta_file)),])
# Load the tab-delimited file containing the values to be analyzed (samples names in the first column)
otu_file <- read.csv(file = input_otu,sep = '\t',row.names = 1,header = T, stringsAsFactors = F, check.names = F)
# Clean table from empty lines
otu_file <- otu_file[!apply(is.na(otu_file) | otu_file =="",1,all),]


############ Express each sample at the selected taxonomy level ##############################
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

############ Convert files to desirable format ######################################

# keep only those rows that appear in the mapping file
otu_file <- otu_file[,rownames(meta_file)]
# OTU-table and mapping file should have the same order and number of sample names
# Order the OTU-table by sample names (ascending)
otu_file <- otu_file[,order(names(otu_file))]
# Transpose OTU-table and convert format to a data frame
otu_file <- data.frame(t(otu_file))
# keep only those that appear in the otu file
meta_file <- meta_file[rownames(otu_file),]

############ Checking for and installing packages ####################################

packages <-c("ade4","dplyr","GUniFrac","phangorn","cluster","fpc","markovchain", 'spgs','caret','nnet','gtools') 
# Function to check whether the package is installed
InsPack <- function(pack)
{
  if ((pack %in% installed.packages()) == FALSE) {
    install.packages(pack,dependencies = T,quiet = T)
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

############ Divide the dataset into different timepoints ########################

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


############ Clustering of samples on all timepoints ###############################

# Initialize a list to save Bayesian Information Criterion for all the timepoints 
BIC_list = list()

# Setting the colour code for the exported plots
xrwmata = c('saddlebrown','cyan3','olivedrab4','sienna1','orange2','yellowgreen','violetred2','rosybrown','orchid4','salmon3')


# Calculate the UniFrac distance matrix for comparing microbial communities
for (name in names(timepoint_list)){
  unifracs <- GUniFrac(otu.tab = timepoint_list[[name]] ,tree = rooted_tree, alpha = c(0.0,0.5,1.0))$unifracs
  
  # Weight on abundant lineages so  the distance is not dominated by highly abundant lineages with 0.5 having the best power
  unifract_dist <- unifracs[, , "d_0.5"]
  
  PAM_clustering <- function(unifract_dist,k){

    clusters = pam(x = unifract_dist, k = k,diss = T)$clustering
    medoids = pam(x = unifract_dist, k = k,diss = T)$medoids
    avg_width = pam(x = unifract_dist, k = k,diss = T)$silinfo$clus.avg.widths

    return (list(clusters,medoids,avg_width))
  }
  
  
  optimal_k<- function(unifract_dist, clustering_method){
    calinski_harabasz_values <-c()
  
    for (k in 2:9){
      clustering_results = PAM_clustering(unifract_dist = unifract_dist, k = k)
      clusteringP <- clustering_results[[1]]
      medoids = clustering_results[[2]]
      avg_width = clustering_results[[3]]
      calinski_harabasz_values= c(calinski_harabasz_values,(calinhara(x = unifract_dist,cn = k, clustering = clusteringP)))
    }
    
    best_delta_score<- which.min(diff(calinski_harabasz_values))
    highest <- which.max(calinski_harabasz_values) +1
    best_final_score <- calinski_harabasz_values[best_delta_score]- calinski_harabasz_values[highest] - min(diff(calinski_harabasz_values))
    if (best_final_score>0){
      best_final_score <- which.max(calinski_harabasz_values) +1
    }
    else {
      best_final_score <- which.min(diff(calinski_harabasz_values)) +1
    }
    return (list(best_final_score, calinski_harabasz_values))
  }
  
  # Set the best k for the timepoint
  optimal_k_results = optimal_k(unifract_dist = unifract_dist,clustering_method = clustering_method)
  best_k = optimal_k_results[[1]]
  calinski_harabasz_values = optimal_k_results[[2]]
  
  # Export plot with the Calinski-Harabasz scores for each k
  jpeg(filename =paste(dir_with_plots, paste(paste('Calinski-Harabasz_index',name,sep = ' of Timepoint '),'jpeg',sep = ' .'),sep = '/'))
  plot (calinski_harabasz_values,type = 'l' ,x = 2:9, col='blue',main = 'Calinski-Harabasz scores of different #Clusters' , ylim=(c(0,max(calinski_harabasz_values) + max(calinski_harabasz_values)*0.2)))
  legend("topright", c("PAM"), fill=c("blue"))
  dev.off()
  
  # Perform clustering
  clustering_results <- PAM_clustering(unifract_dist = unifract_dist ,k =best_k)
  clusters = clustering_results[[1]]
  medoids = clustering_results[[2]]
  avg_width = clustering_results[[3]]
  # Assing the samples to the estimated clusters
  samples_on_clusters[meta_file[row.names(unifract_dist),'Sample'],name]<- clusters
  samples_on_clusters[samples_on_clusters==0] <- NA
  
  # Export MDS plot of the samples on the timepoint
  scall <- cmdscale(d = unifract_dist,eig = F,k = 2)
  jpeg(filename = paste(dir_with_plots,paste('MDS Plot of',paste('Timepoint',name,sep=' '),sep = ' '),sep = '/'))
  plot(scall , main = name, col = xrwmata[clusters] , pch = clusters)
  dev.off()
  
  #############Gaussian Mixture Model Based Test to evaluate whether the dataset forms groups ###############################################
  
  gmm_testing = Mclust(data = as.dist(unifract_dist),G = c(1,max(samples_on_clusters[,name], na.rm = T)),  modelNames = c("EII","VII","EEI","EVI","VEI","VVI") , verbose = F)
  
  BIC_list[[name]] = gmm_testing$BIC
  if (gmm_testing$G == 1){
    samples_on_clusters[,name] = rep(1,nrow(samples_on_clusters))
  }
  
}
