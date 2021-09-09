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
dir_with_plots= 'Calinski-Harabasz Plots'     # <---- CHANGE ACCORDINGLY
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
namesrow <- c()
namesample <-c()
for (i in rownames(meta_file)){
  namesrow <- append(namesrow, paste(ifelse(nchar(i)==5,yes = 'X',no = 'X0'), i, sep = ""))
  namesample <- append(namesample, paste(ifelse(nchar(meta_file[i,'Sample'])==2,yes = "X0",no = "X"),meta_file[i,'Sample'], sep = "")) 
}
rownames(meta_file) <- namesrow
meta_file[,'Sample'] <- namesample
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

packages <-c("ade4","dplyr","GUniFrac","phangorn","cluster","fpc","markovchain", 'spgs') 
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
################ CLUSTERING SAMPLES ON TIMEPOINTS ########################################
##########################################################################################
# Calculate the UniFrac distance matrix for comparing microbial communities
for (name in names(timepoint_list)){

unifracs <- GUniFrac(otu.tab = timepoint_list[[name]] ,tree = rooted_tree, alpha = c(0.0,0.5,1.0))$unifracs

# Weight on abundant lineages so the distance is not dominated by highly abundant lineages with 0.5 having the best power
unifract_dist <- unifracs[, , "d_0.5"]

PAM_clustering <- function(unifract_dist,k){
  return (pam(x = unifract_dist, k = 4,diss = T)$clustering)
}

H_clustering <- function(unifract_dist,k){
  return(cutree(hclust(d = as.dist(unifract_dist),method = 'ward.D2'),k=k))
}

optimal_k<- function(unifract_dist, clustering_method){
  calinski_harabasz_values <-c()
  if (clustering_method == 'Hierarchical'){
    for (k in 2:9){
      clustering <- H_clustering(unifract_dist = unifract_dist, k = k)
      calinski_harabasz_values= c(calinski_harabasz_values,(calinhara(x = unifract_dist,cn = k, clustering = clustering)))
    }
  }
  else if (clustering_method == 'PAM'){
  for (k in 2:9){
    clustering <- PAM_clustering(unifract_dist = unifract_dist, k = k)
    calinski_harabasz_values= c(calinski_harabasz_values,(calinhara(x = unifract_dist,cn = k, clustering = clustering)))
    }
  }
  kaliterotero<- which.min(diff(calinski_harabasz_values))+1
  highest <- which.max(calinski_harabasz_values)
  kalitero <- calinski_harabasz_values[kaliterotero]- calinski_harabasz_values[highest] - min(diff(calinski_harabasz_values))
  if (kalitero>0){
    kalitero <- which.max(calinski_harabasz_values) +1
  }
  else {
    kalitero <- which.min(diff(calinski_harabasz_values))+1
  }
  dir.create(dir_with_plots, showWarnings = F)
  jpeg(filename =paste(dir_with_plots, paste(paste('Calinski-Harabasz_index',name,sep = ' of Timepoint '),'jpeg',sep = ' .'),sep = '/'))
  plot(calinski_harabasz_values,x = 2:9)
  dev.off()
  return (kalitero)
}

k=optimal_k(unifract_dist = unifract_dist,clustering_method = clustering_method)


if (clustering_method == 'Hierarchical'){
  clusters <- H_clustering(unifract_dist = unifract_dist, k = k )
}
else if (clustering_method == 'PAM'){
  clusters <- PAM_clustering(unifract_dist = unifract_dist ,k =k)
}
samples_on_clusters[meta_file[row.names(unifract_dist),'Sample'],name]<- clusters
samples_on_clusters[samples_on_clusters==0] <- NA
}

######################### END OF SECTION #################################################

##########################################################################################
######################## SPECIFY TAXA ON CLUSTERS   ######################################
##########################################################################################


taxa_per_cluster <- function(taxa_matrix,samples_on_clusters,timepoint_list){
  taxa_clusters <- list()
  for (i in names(timepoint_list)){
    for (j in 1:max(samples_on_clusters[,i],na.rm = T)){
      taxa_clusters[[paste(as.character(i),as.character(j),sep = ' cluster ')]] <- as.matrix(x = apply(X = taxa_matrix[paste(rownames(samples_on_clusters[samples_on_clusters[,i]==j,]),i,sep = ''),], MARGIN = 2, FUN = representation_method),decreasing = T)      
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
################ WRITE TAB DELIMITED FILES WITH THE OUTPUTS  #############################
##########################################################################################


dir.create(dir_with_files, showWarnings = F)


colnames(taxa_clusters)<-lapply(X = colnames(taxa_clusters), FUN = function(x){paste('Timepoint',x,sep = ' ')})
row.names(samples_on_clusters) <- lapply(X = rownames(samples_on_clusters), FUN = function(x){substr(x,start=2,stop= 4)})

write.csv(x = samples_on_clusters, file = paste(dir_with_files, "Samples_in_Timepoint-specific_Clusters.csv",sep = '/'), row.names = T)
write.csv(x = taxa_clusters,       file = paste(dir_with_files, "Taxonomic_profile_of_clusters.csv", sep = '/'), row.names = T)

######################### END OF SECTION #################################################

################### MARKOVIAN CHAIN CHECK ################################################
state_calculator = function(infants,t0,t2){
  states <- c()
  for (j in 1:nrow(unique(infants[,t0:t2]))){
    if (!any(is.na(unique(infants[,t0:t2])[j,]))){
    states[[j]] <- as.vector (unique(infants[,t0:t2])[j,])
    }
  }
  return (states)
}  

counting_states = function (infants,t0,t2,states) {
  counter <- c(rep(0,length(states)))
  for (i in 1:nrow(infants)){
    if (!any(is.na(unique(infants[i,t0:t2])))){
      for (j in 1:length(states)){
        if (all(infants[i,t0:t2]==states[[j]])){
          counter[j]= counter[j]+1
          } 
      }
    }
  }
  return (counter)
}

Ss <- c()
dof <- c()
for (t0 in 1:(ncol(infants)-2)){
  t1 <- t0+1
  t2 <- t0+2
  ijk_states <- state_calculator(infants = infants, t0 = t0 , t2 = t2 )
  ijk_counter <- counting_states(infants = infants ,t0 = t0, t2 = t2, states = ijk_states)
  
  jk_states <-lapply(ijk_states, function(x){x[2:3]})
  jk_counter <- counting_states(infants = infants ,t0 = t1,t2 = t2, states = jk_states)
  ij_states <- lapply(ijk_states, function(x){x[1:2]})
  
  p_jk <- jk_counter/sum(jk_counter)
  S =0
  for (i in 1:length(jk_counter)){
      S = S + (ijk_counter[i]-jk_counter[i]*p_jk[[i]])**2/jk_counter[i]*p_jk[[i]]
  }
  degrees_of_freedom = length(unique(infants[,t0])) - length(unique(ij_states)) + length(unique(ijk_states)) +1  
  Ss[t0] <- S
  dof [t0] <- degrees_of_freedom
}
Markovian_Property_Per_Timepoints <-c()
for (i in 1: length(dof)){
  Markovian_Property_Per_Timepoints[i] = Ss[i] < qchisq(.95,df = dof[i])
}
Markovian_Property = all(Markovian_Property_Per_Timepoints)
Markovian_Property

################### MARKOVIAN CHAIN CHECK 2 ##############################################
infants<- samples_on_clusters[is.na(samples_on_clusters[,adult_timepoint_name]),1:ncol(samples_on_clusters)-1]

independed_clusters <- cumsum(as.vector(apply(X = infants,MARGIN = 2,FUN = function(x){max(x,na.rm = T)})))
for (i in 2:ncol(infants)){
  for (j in 1:max(infants[,i],na.rm = T)){
     infants[infants[,i]==j,i] <- independed_clusters[i-1]+j
  }
}
apply(samples_on_clusters[is.na(samples_on_clusters[,adult_timepoint_name]),1:ncol(samples_on_clusters)-1], 2, function(x){max(x,na.rm = T)})
apply(infants,MARGIN = 2,FUN = function(x){max(x,na.rm = T)})

testaki <- markov.test(infants)
testaki$statistics["X-squared"]




##########################################################################################
################ PLOT PRINCIPAL COMPONENTS FOR ALL TIMEPOINTS ############################
##########################################################################################

normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}
for (tt1 in names(timepoint_list)){
  dir.create('Principal Component Analysis',showWarnings = F)
  tt2= timepoint_list[[tt1]][,colSums(timepoint_list[[tt1]])!=0]
  tt2 = normalize(tt2)
  tt = prcomp(x = tt2,scale. = F, center = T)
  ttt = summary(tt)
  jpeg(filename = paste('Principal Component Analysis',paste(paste('Timepoint', tt1, sep = ' '),'jpeg',sep = '.'),sep = '/'))
  plot(ttt$x[,1],ttt$x[,2],main = paste('Timepoint', tt1, sep = ' ') ,xlab = paste ((ttt$importance[2,1]*100),'%', sep = ' '), ylab = paste((ttt$importance[2,2]*100),'%',sep = ' '))
  dev.off()
}

######################### END OF SECTION #################################################


############## PRACTICE ############################



#lapply(rownames(samples_on_clusters),FUN = function(x){paste('X',x, sep = "")})


#meta_file[meta_file[,'Sample']== unlist(lapply(rownames(samples_on_clusters),FUN = function(x){paste('X',x, sep = "")}))]
#otu_file[paste(meta_file[meta_file[,'Timepoint']=='01','Sample'], '01', sep = ""),1:5]
#samples_on_clusters[substr(x = meta_file[meta_file[,'Timepoint']=='01','Sample'], start = 2,stop = 5),1]

apply(samples_on_clusters, 2, FUN = function(x){max(x,na.rm = T)})
############## COMMENTS ############################