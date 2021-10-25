##########################################################################################
############ ON THIS SECTION YOU CAN SELECT YOUR DESIRED PARAMETERS ######################
##########################################################################################

########### PLEASE FOLLOW THE INSTRUCTIONS CAREFULLY #####################################

# Please set the directory of the script as the working folder (e.g D:/studyname/Data/Chronos/)
# Note: the path is denoted by forward slash "/".
working_directory = "~/Working_Chronos/Cronos_Final/"  #<--- CHANGE ACCORDINGLY !!!

# Here we set the working directory as you selected in order to find your files and R scripts
# Please do not interfere with the next line

setwd(working_directory)


# Please give the file name of the normalized OTU-table without taxonomic classification
input_otu = "SOTUs-Table.tab"           #<--- CHANGE ACCORDINGLY !!!
# Please give the name of the meta-file that contains individual sample information
input_meta = "Mapping_File_Inf_St.csv"                #<--- CHANGE ACCORDINGLY !!!
# Please give the name of the phylogenetic tree constructed from the OTU sequences
input_tree = "SOTUs-NJTree.tre"         #<--- CHANGE ACCORDINGLY !!!


# Please specify if the file contains both adult and infant data:
# If it contains adult data specify the name of the column (i.e. timepoint name) where they are saved
# If it doesn't, leave it blank (like this: '' ).
adult_timepoint_name = 'MM'              #<--- CHANGE ACCORDINGLY


# Please select the taxon in which the samples will be analyzed
# Either type it in ' ' e.g. 'Order' or select a number between 1 and 5, where:
# 1: Domain,
# 2: Phylum
# 3: Class
# 4: Order
# 5: Family

taxonomic_level='Family'                          # <---- CHANGE ACCORDINGLY


# Please specify the number of iterations to use for the stratified test split on modeling.
# Final Stratified accuracy of the modeling will derive as the mean accuracy from
# the number selected. 
# This number is strongly linked to the time Chronos will need to run.
# Greater numbers lead to more time needed (default is 20)
splitting_times = 2                           # <---- CHANGE ACCORDINGLY

# Please select the action Cronos should do if the same parameters are already used for analysis before
# It can be either 'Continue' or 'Stop'.
action = 'Stop'                           # <---- CHANGE ACCORDINGLY


########### Now that you selected the parameters you can select everything (Ctrl+A) and press run (Ctrl+Enter)  ##########

##########################################################################################################################
###################### PLEASE DO NOT CHANGE ANYTHING BEYOND THIS POINT ###################################################
##########################################################################################################################

##########################################################################################################################
#################################### CREATING OUTPUT DIRECTORY ###########################################################
##########################################################################################################################

# Setting the name of the directory where Cronos outputs are stored 
output_dir = paste('Cronos',date(),sep =' ')
# Create a log file with all the parameters used
parameters = matrix('',ncol = 1 ,nrow = 6)
rownames(parameters) = c('input_meta','input_tree','input_otu','taxonomic_level','adult_timepoint_name','splitting_times')
parameters[,1]= c(input_meta,input_tree,input_otu,taxonomic_level,adult_timepoint_name,splitting_times)

# Checking whether Cronos was run before with the same parameters
new_run = T
for (directory in list.dirs()[2:length(list.dirs())]){
  previous_runs = list.files(path = directory)
  existing_runs = read.csv(paste(directory,previous_runs[grepl(pattern = '_log.csv', x = previous_runs , ignore.case = F)], sep = '/'))
  if (all(existing_runs[,2] == parameters[,1])){
    new_run = F
  }
}


if (new_run==T || (new_run == F & action =='Continue')){
  
  # Create the directory where Cronos outputs will be stored
  dir.create(output_dir,showWarnings = F)
  # Write file with the parameters on this run
  write.csv(x = parameters,file = paste(output_dir,'Cronos_log.csv',sep = '/') ,row.names = T, na = ' ')
  
##################################### SECTION ############################################################################
#################################### CLUSTERING ##########################################################################
##########################################################################################################################

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

packages <-c("ade4","dplyr","GUniFrac","phangorn","cluster","fpc","markovchain", 'spgs','caret','nnet','gtools', 'mclust','igraph') 
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
colours_ploting = c('saddlebrown','cyan3','olivedrab4','sienna1','orange2','yellowgreen','violetred2','rosybrown','orchid4','salmon3', colors())

# Initialize a list to save the medoid information
medoids = matrix(NA, nrow = 9, ncol = length(timepoint_list))
colnames(medoids) = names(timepoint_list)

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
  jpeg(filename =paste(output_dir, paste(paste('Calinski-Harabasz_index',name,sep = ' of Timepoint '),'jpeg',sep = ' .'),sep = '/'))
  plot (calinski_harabasz_values,type = 'l' ,x = 2:9, col='blue',main = 'Calinski-Harabasz scores of different #Clusters' , ylim=(c(0,max(calinski_harabasz_values) + max(calinski_harabasz_values)*0.2)))
  legend("topright", c("PAM"), fill=c("blue"))
  dev.off()
  
  # Perform clustering
  clustering_results <- PAM_clustering(unifract_dist = unifract_dist ,k =best_k)
  clusters = clustering_results[[1]]
  medoids_temp = rep(NA,9-length(clustering_results[[2]]))
  medoids[,name] = c(clustering_results[[2]], medoids_temp)
  #  medoids = c(medoids, clustering_results[[2]])
  avg_width = clustering_results[[3]]
  # Assing the samples to the estimated clusters
  samples_on_clusters[meta_file[row.names(unifract_dist),'Sample'],name]<- clusters
  samples_on_clusters[samples_on_clusters==0] <- NA
  
  # Export MDS plot of the samples on the timepoint
  scall <- cmdscale(d = unifract_dist,eig = F,k = 2)
  jpeg(filename = paste(output_dir,paste('MDS Plot of',paste('Timepoint',name,sep=' '),sep = ' '),sep = '/'))
  plot(scall , main = name, col = colours_ploting[clusters] , pch = clusters)
  dev.off()
  
  #############Gaussian Mixture Model Based Test to evaluate whether the dataset forms groups ###############################################
  
  gmm_testing = Mclust(data = as.dist(unifract_dist),G = c(1,max(samples_on_clusters[,name], na.rm = T)),  modelNames = c("EII","VII","EEI","EVI","VEI","VVI") , verbose = F)
  
  BIC_list[[name]] = gmm_testing$BIC
  if (gmm_testing$G == 1){
    samples_on_clusters[,name] = rep(1,nrow(samples_on_clusters))
  }
  
}


######################################### SECTION ########################################################################
#################################### TRANSITION ANALYSIS #################################################################
##########################################################################################################################

############ Markovian property check ################################################

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

if (nchar(adult_timepoint_name)>0){
  tps = as.character(sort(as.numeric(colnames(samples_on_clusters)[colnames(samples_on_clusters)!='ot' & colnames(samples_on_clusters)!=adult_timepoint_name]),decreasing = F))
} else {
  tps = as.character(sort(as.numeric(colnames(samples_on_clusters[colnames(samples_on_clusters)!='ot'])),decreasing = F))
}


infants = samples_on_clusters[,tps]
infants = infants[!apply(infants,1,function(x){all(is.na(x))}),]
infants_on_clusters <- infants
independed_clusters <- cumsum(as.vector(apply(X = infants,MARGIN = 2,FUN = function(x){max(x,na.rm = T)})))

for (i in 2:ncol(infants_on_clusters)){
  for (j in 1:max(infants_on_clusters[,i],na.rm = T)){
    infants_on_clusters[infants_on_clusters[,i]==j,i] <- independed_clusters[i-1]+j
  }
}

Markov_Test_X2 <- function(infants){
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
    ij_counter <- counting_states(infants = infants ,t0 = t0, t2 = t1,states = ij_states)
    
    p_jk <- jk_counter/sum(jk_counter)
    S =0
    for (i in 1:length(jk_counter)){
      S = S + (ijk_counter[i]-ij_counter[i]*p_jk[[i]])**2/ij_counter[i]*p_jk[[i]]
    }
    degrees_of_freedom = length(unique(infants[,t0])) - length(unique(ij_states)) + length(unique(ijk_states)) +1  
    Ss[t0] <- S
    dof [t0] <- degrees_of_freedom
  }
  Markovian_Property_Per_Timepoints <-c()
  for (i in 1: length(dof)){
    Markovian_Property_Per_Timepoints[i] = Ss[i] > qchisq(.05,df = dof[i],lower.tail = F)
  }
  Markovian_Property_of_transitions <- matrix(NA, nrow = 3,ncol = length(dof))
  rownames(Markovian_Property_of_transitions) = c('Markov_Property_per_transition','X2', 'Degrees_of_freedom')
  Markovian_Property_of_transitions [1,] <- Markovian_Property_Per_Timepoints
  Markovian_Property_of_transitions [2,] <- Ss
  Markovian_Property_of_transitions [3,] <- dof
  return (Markovian_Property_of_transitions)
}

Markovian_Check <- Markov_Test_X2(infants = infants)
Markovian_property= all(Markovian_Check[1,]==1)

############ Claim transition matrix ################################################

matrix_of_transitions <- matrix(0, nrow = max(infants_on_clusters, na.rm = T), ncol = max(infants_on_clusters,na.rm = T))
for (i in 1:(ncol(infants_on_clusters)-1)){
  for (j in min(infants_on_clusters[,i],na.rm = T):max(infants_on_clusters[,i],na.rm = T)){
    for (k in min(infants_on_clusters[,i+1],na.rm = T):max(infants_on_clusters[,i+1],na.rm = T)){
      
      matrix_of_transitions[j,k] = sum(infants_on_clusters[infants_on_clusters[,i]==j,(i+1)]==k,na.rm = T)/sum(infants_on_clusters[,i]==j & !is.na(infants_on_clusters[,i+1]),na.rm = T)
      
    }
  }
}

############ Specify taxa on clusters #################################################

taxa_per_cluster <- function(taxa_matrix,samples_on_clusters,timepoint_list, representation_method){
  taxa_clusters <- list()
  for (i in colnames(samples_on_clusters)){
    for (j in min(samples_on_clusters[,i], na.rm = T):max(samples_on_clusters[,i],na.rm = T)){
      timiclj = samples_on_clusters[!is.na(samples_on_clusters[,i]),i]
      metatimepoint= meta_file[meta_file[,"Timepoint"]==i,]
      
      taxa_clusters[[paste(as.character(i),as.character(j),sep = ' cluster ')]]= taxa_matrix[medoids[j,i],]
      
    }
  }
  clustering_taxa <- matrix(0,  nrow = ncol(taxa_matrix), ncol = length(taxa_clusters))
  row.names(clustering_taxa)<- colnames(taxa_matrix)
  colnames(clustering_taxa)<- names(taxa_clusters)
  
  for (i in colnames(clustering_taxa)){
    for (j in rownames(clustering_taxa)){
      clustering_taxa[j,i] <- taxa_clusters[[i]][,j]
    }
  }
  return (clustering_taxa)
}

taxa_clusters <- taxa_per_cluster(taxa_matrix = taxa_matrix, samples_on_clusters = samples_on_clusters, timepoint_list = timepoint_list, representation_method = representation_method)

samples_on_clusters = samples_on_clusters[,names(timepoint_list)[names(timepoint_list)!='ot']]
samples_on_clusters = samples_on_clusters[,c(as.character(sort(as.numeric(colnames(samples_on_clusters)[colnames(samples_on_clusters)!=adult_timepoint_name]),decreasing = F,na.last = T)),adult_timepoint_name)]


infants<- samples_on_clusters[is.na(samples_on_clusters[,adult_timepoint_name]),1:ncol(samples_on_clusters)-1]

############ Write comma delimited files with outputs   ####################################

colnames(taxa_clusters)<-lapply(X = colnames(taxa_clusters), FUN = function(x){paste('Timepoint',x,sep = ' ')})

write.csv(x = samples_on_clusters, file = paste(output_dir, "Samples_in_Timepoint-specific_Clusters.csv",sep = '/'), row.names = T)
write.csv(x = taxa_clusters,       file = paste(output_dir, "Taxonomic_profile_of_clusters.csv", sep = '/'), row.names = T)


No_clusters_per_timepoint = apply(infants, 2,function(x){max(x,na.rm = T)})
transition_names = c()
for (j in 1:length(No_clusters_per_timepoint)){
  for (i in 1:(No_clusters_per_timepoint[j])){
    transition_names= c(transition_names,paste( paste('Timepoint',colnames(infants)[j],sep = ' '), paste( 'Cluster',i,sep = ' '),sep = ' '))
  }
}
rownames(matrix_of_transitions)= paste('From ', transition_names,sep = '')
write.csv(x = matrix_of_transitions,file = paste(output_dir,'Transition_Matrix.csv',sep = '/') , row.names = T)



############ Plotting the transitions ####################################

transition_graph <- graph_from_adjacency_matrix(adjmatrix = matrix_of_transitions, weighted = T ,mode = 'directed' )

graph.colours = c()
for (i in 1:ncol(infants)){
  graph.colours = c(graph.colours,rep(i,No_clusters_per_timepoint[i]))
}
jpeg(filename = paste(output_dir,'Transition_graph',sep = '/'))
plot(x = transition_graph, vertex.color = colours_ploting[graph.colours], vertex.shape = 'circle', layout= layout_as_tree ,edge.arrow.size = transition_graph[])
legend("topleft", paste('Clusters at Timepoint',colnames(infants),sep = ' '), fill = colours_ploting[1:ncol(infants)])
dev.off()


##################################### SECTION ############################################################################
#################################### MODELING ############################################################################
##########################################################################################################################


###################### Initializing Lists to save metrics of modeling ############################
# Initializing lists to write accuracies of models
Train_sets_LOO <- c()
Test_sets_LOO  <- c()
Train_sets_stratified_split <- c()
Test_sets_stratified_split  <- c()

# Naming the effectors to be imputed on the model
effectors = colnames(meta_file)[3:ncol(meta_file)]

###################### Perform Tasks on all Timepoints ###############################

# Loop to create a model for each timepoint on the dataset
for (end_timepoint in head(as.numeric(rev(colnames(infants))),n = (ncol(infants))-1)){
  
  
  # Setting the timepoints from which the different models derive  
  timepoints_to_perform = rev(colnames(infants)[as.numeric(colnames(infants)) < end_timepoint])
  
  ###################### Leave one out Null #######################################################
  
  # Calculate mean accuracy of prediction on final timepoint, via Leave One Out method of splitting Train-Test sets
  LOO_multilogreg_Null <- function (infants,timepoint,end_timepoint,meta_file){
    
    logreg <- as.data.frame(matrix(NA,nrow = nrow(infants),ncol = 2))
    colnames(logreg)= c(paste('Cluster_at',end_timepoint,sep = '_'), timepoint)
    rownames(logreg)= rownames(infants)
    
    for (modeling_sample in (rownames(infants))){
      
      logreg[modeling_sample, 2] = ifelse(test = length(infants[modeling_sample,timepoint]), yes = infants[modeling_sample,timepoint], no = NA)
      logreg[modeling_sample, 1] = ifelse(test = length(infants[modeling_sample,as.character(end_timepoint)]), yes = infants[modeling_sample,as.character(end_timepoint)], no = NA)
      
    }
    
    logreg = logreg[complete.cases(logreg),]
    
    acc <- c()
    tra <- c()
    form = as.formula(paste(paste('Cluster_at',end_timepoint,sep = '_'),'~.',sep = ''))
    for (j in 1:nrow(logreg)){
      trainset = logreg[1:nrow(logreg)!=j,]
      testset = logreg[j,]
      
      prediction_model <- multinom(formula = form , data = trainset ,censored = T, model = T)
      
      training_set_accuracies <- prediction_model %>% predict(trainset)
      test_set_predictions <- prediction_model %>% predict(testset)
      acc[j] = (sum(test_set_predictions == testset[,1])/ nrow(testset)) *100
      tra[j] = (sum(test_set_predictions == trainset[,1])/nrow(trainset))*100
      
    }
    
    return (c(mean(tra),mean(acc)))
  }
  
  ###################### Stratified Train/Test splits Null ########################################
  
  # Calculate accuracy of prediction on final timepoint, as a mean of 100 accuracies with different
  # Train-Test stratified splits. Stratfication is performed on all metadata categories 
  multilogreg_stratified_Null <- function(infants,optimal_splitting_percentage,timepoint,end_timepoint,meta_file, splitting_times){
    
    logreg <- as.data.frame(matrix(NA,nrow = nrow(infants),ncol = 2))
    colnames(logreg)= c(paste('Cluster_at',end_timepoint,sep = '_'),timepoint)
    rownames(logreg)= rownames(infants)
    
    for (modeling_sample in (rownames(infants))){
      
      logreg[modeling_sample, 2] = ifelse(test = length(infants[modeling_sample,as.character(timepoint)]), yes = infants[modeling_sample,as.character(timepoint)], no = NA)
      logreg[modeling_sample, 1] = ifelse(test = length(infants[modeling_sample,as.character(end_timepoint)]), yes = infants[modeling_sample,as.character(end_timepoint)], no = NA)
      
    }
    logreg = logreg[complete.cases(logreg),]
    
    train_index <- createDataPartition(y = logreg[,paste('Cluster_at',end_timepoint,sep = '_')], p = optimal_splitting_percentage, list = F, times = splitting_times, groups = max(apply(X = logreg, MARGIN = 2, FUN = function(x){length(unique(x))})) )
    
    acc <- c()
    tra <- c()
    weighted_modeling_factors <- list()
    fo = as.formula(paste(paste('Cluster_at',end_timepoint,sep = '_'),'~.',sep = ''))
    for (j in 1:ncol(train_index)){
      trainset = logreg[train_index[,j],]
      testset = logreg[-train_index[,j],]
      if (all(apply(X = logreg, MARGIN = 2, FUN = function(x){length(unique(x))})== apply(X = trainset, MARGIN = 2, FUN = function(x){length(unique(x))}))){
        
        prediction_model <- multinom(formula = fo, data = trainset ,censored = T, model = T)
        
        training_set_accuracies <- prediction_model %>% predict(trainset)
        test_set_predictions <- prediction_model %>% predict(testset)
        acc[j] = (sum(test_set_predictions == testset[,1])/ nrow(testset)) *100
        tra[j] = (sum(training_set_accuracies == trainset[,1])/nrow(trainset))*100
      }
    }
    
    acc = na.omit(acc)
    tra = na.omit(tra)
    
    Test_set_stratified_accuracy <- mean(acc)
    Train_set_stratified_accuracy = mean(tra)
    trainsd <- sd(tra)
    testsd <- sd(acc)
    
    return (c(Test_set_stratified_accuracy,Train_set_stratified_accuracy, trainsd,testsd))
  }
  
  ###################### Calculate optimal split percentage ############################################
  # Calculate the optimal split percentage between .65 and .90 for predicting the final timepoint on the
  # exact previous one to be used in all splits for all timepoints. The criterion is the best accuracy  
  # as derived by the previous function with argument the different percentages.
  optimal_splitting_percentage <- function (infants){
    all_calculated_accuracies <- c()
    best_splitting_percentage <- c()
    for (percentage_to_calculate in 65:90){
      percentage_to_calculate = percentage_to_calculate/100
      telika <- multilogreg_stratified_Null(infants = infants, optimal_splitting_percentage = percentage_to_calculate, timepoint = tail(timepoints_to_perform,1), end_timepoint = end_timepoint,meta_file = files,splitting_times = splitting_times)
      Test_set_stratified_accuracy = telika[1]
      Train_set_stratified_accuracy = telika[2]
      trainsd = telika[3]
      testsd = telika[4]
      all_calculated_accuracies= c(all_calculated_accuracies,Test_set_stratified_accuracy)
      best_splitting_percentage = c(best_splitting_percentage,percentage_to_calculate)
    }
    
    return (best_splitting_percentage[which.max(all_calculated_accuracies)])
  }
  
  # Assign the best percentage on a variable to use on the final calculations
  optimal_splitting_percentage <- optimal_splitting_percentage(infants = infants)
  
  ###################### Initializing variables Null ########################################
  Test_set_LOO_accuracy <- c()
  Train_set_LOO_accuracy <-c()
  Test_set_stratified_accuracy <- c()
  Train_set_stratified_accuracy <- c()
  trainsd <- c()
  testsd <- c()
  weighted_modeling_factors <- list()
  
  ###################### Calculate Accuracies on timepoints Null ############################
  
  # Loop to create a model from every timepoint to the final timepoint as selected from the first loop
  for (timepoint in timepoints_to_perform){
    
    LOOTeliko <- LOO_multilogreg_Null(infants = infants, timepoint = timepoint,end_timepoint = end_timepoint,meta_file = files)
    
    Train_set_LOO_accuracy = c(Train_set_LOO_accuracy,LOOTeliko[1])
    Test_set_LOO_accuracy = c(Test_set_LOO_accuracy,LOOTeliko[2])
    
    MultiLogReg_results <- multilogreg_stratified_Null(infants = infants, optimal_splitting_percentage = optimal_splitting_percentage,timepoint =  timepoint,meta_file =  files, end_timepoint = end_timepoint,splitting_times = splitting_times)
    
    
    Test_set_stratified_accuracy <- c(Test_set_stratified_accuracy,MultiLogReg_results[1])
    Train_set_stratified_accuracy <- c(Train_set_stratified_accuracy, MultiLogReg_results[2])
    trainsd <- c(trainsd, MultiLogReg_results[3])
    testsd <- c(testsd, MultiLogReg_results[4])
  }  
  
  
  ###################### Save accuracies on the lists ######################################
  
  # Adding to the lists in order to export the correct matrix
  Train_sets_LOO = rbind (Train_sets_LOO, c('Null',rep(0,(ncol(infants)-length(Train_set_LOO_accuracy)-1)),Train_set_LOO_accuracy))
  Test_sets_LOO  = rbind (Test_sets_LOO,  c('Null',rep(0,(ncol(infants)-length(Test_set_LOO_accuracy) -1)), Test_set_LOO_accuracy))
  Test_sets_stratified_split  = rbind (Test_sets_stratified_split,  c('Null',rep(0,(ncol(infants)-length(Test_set_stratified_accuracy) -1)), Test_set_stratified_accuracy))
  Train_sets_stratified_split = rbind (Train_sets_stratified_split, c('Null',rep(0,(ncol(infants)-length(Train_set_stratified_accuracy)-1)), Train_set_stratified_accuracy))
  
  
  ###################### Perform all Tasks with different combinations of metadata ############
  
  for (Ncomb in 1: (ncol(meta_file)-3)){
    # Create a vector of all possible feature combinations
    effector_combinations = combinations(v = effectors,r = Ncomb, repeats.allowed = F, n = length(effectors))
    # Loop to select all possible combinations of features from the dataset  
    for (combination in 1:nrow(effector_combinations)){
      # Select the table with the corresponding features on the corresponding timepoint
      files = meta_file[,c('Sample','Timepoint',effector_combinations[combination,])]
      
      ###################### Leave one out Combinations #######################################################
      
      # Calculate mean accuracy of prediction on final timepoint, via Leave One Out method of splitting Train-Test sets
      LOO_multilogreg_One <- function (infants,timepoint,end_timepoint,meta_file,Ncomb){
        
        logreg <- as.data.frame(matrix(NA,nrow = nrow(infants),ncol = ncol(meta_file)))
        colnames(logreg)= c(paste('Cluster_at',end_timepoint,sep = '_'), timepoint , colnames(meta_file)[3:ncol(meta_file)])
        rownames(logreg)= rownames(infants)
        
        for (modeling_sample in (rownames(infants))){
          subTab = meta_file[meta_file[,1] %in% modeling_sample & meta_file[,2] == as.numeric(end_timepoint), 3:ncol(meta_file)]
          if (Ncomb > 1 ){
            if(nrow(subTab) == 1){
              logreg[modeling_sample, colnames(meta_file)[3:ncol(meta_file)]] = subTab
            }
          }
          else if (Ncomb == 1){
            if (length(subTab) == 1){
              logreg[modeling_sample, colnames(meta_file)[3:ncol(meta_file)]] = subTab
            }
          }
          logreg[modeling_sample, 2] = ifelse(test = length(infants[modeling_sample,timepoint]), yes = infants[modeling_sample,timepoint], no = NA)
          logreg[modeling_sample, 1] = ifelse(test = length(infants[modeling_sample,as.character(end_timepoint)]), yes = infants[modeling_sample,as.character(end_timepoint)], no = NA)
        }
        
        logreg = logreg[complete.cases(logreg),]
        
        acc <- c()
        tra <- c()
        form = as.formula(paste(paste('Cluster_at',end_timepoint,sep = '_'),'~.',sep = ''))
        for (j in 1:nrow(logreg)){
          trainset = logreg[1:nrow(logreg)!=j,]
          testset = logreg[j,]
          
          prediction_model <- multinom(formula = form , data = trainset ,censored = T, model = T)
          
          training_set_accuracies <- prediction_model %>% predict(trainset[,2:ncol(trainset)])
          test_set_predictions <- prediction_model %>% predict(testset[,2:ncol(testset)])
          acc[j] = (sum(test_set_predictions == testset[,1])/ nrow(testset)) *100
          tra[j] = (sum(test_set_predictions == trainset[,1])/nrow(trainset))*100
          
        }
        
        return (c(mean(tra),mean(acc)))
      }
      
      ###################### Stratified Train/Test splits Combinations ########################################
      
      # Calculate accuracy of prediction on final timepoint, as a mean of 100 accuracies with different
      # Train-Test stratified splits. Stratfication is performed on all metadata categories 
      multilogreg_stratified_One <- function(infants,optimal_splitting_percentage,timepoint,end_timepoint,meta_file,Ncomb,splitting_times){
        
        logreg <- as.data.frame(matrix(NA,nrow = nrow(infants),ncol = ncol(meta_file)))
        colnames(logreg)= c(paste('Cluster_at',end_timepoint,sep = '_'),timepoint,colnames(meta_file)[3:ncol(meta_file)])
        rownames(logreg)= rownames(infants)
        
        for (modeling_sample in (rownames(infants))){
          
          subTab = meta_file[meta_file[,1] %in% modeling_sample & meta_file[,2] == as.numeric(end_timepoint), 3:ncol(meta_file)]
          
          if (Ncomb > 1){
            if (nrow(subTab) == 1){  
              logreg[modeling_sample, colnames(meta_file)[3:ncol(meta_file)]] = subTab
            }
          }
          else if (Ncomb == 1){
            if (length(subTab) == 1){
              logreg[modeling_sample, colnames(meta_file)[3:ncol(meta_file)]] = subTab
            }
          }
          logreg[modeling_sample, 2] = ifelse(test = length(infants[modeling_sample,timepoint]), yes = infants[modeling_sample,timepoint], no = NA)
          logreg[modeling_sample, 1] = ifelse(test = length(infants[modeling_sample,as.character(end_timepoint)]), yes = infants[modeling_sample,as.character(end_timepoint)], no = NA)
          
          
        }
        logreg = logreg[complete.cases(logreg),]
        
        train_index <- createDataPartition(y = logreg[,paste('Cluster_at',end_timepoint,sep = '_')], p = optimal_splitting_percentage, list = F, times = splitting_times, groups = max(apply(X = logreg, MARGIN = 2, FUN = function(x){length(unique(x))})) )
        
        acc <- c()
        tra <- c()
        weighted_modeling_factors <- list()
        fo = as.formula(paste(paste('Cluster_at',end_timepoint,sep = '_'),'~.',sep = ''))
        for (j in 1:ncol(train_index)){
          trainset = logreg[train_index[,j],]
          testset = logreg[-train_index[,j],]
          if (all(apply(X = logreg, MARGIN = 2, FUN = function(x){length(unique(x))})== apply(X = trainset, MARGIN = 2, FUN = function(x){length(unique(x))}))){
            
            prediction_model <- multinom(formula = fo, data = trainset ,censored = T, model = T)
            
            training_set_accuracies <- prediction_model %>% predict(trainset[,2:ncol(trainset)])
            test_set_predictions <- prediction_model %>% predict(testset[,2:ncol(testset)])
            acc[j] = (sum(test_set_predictions == testset[,1])/ nrow(testset)) *100
            tra[j] = (sum(training_set_accuracies == trainset[,1])/nrow(trainset))*100
          }
        }
        
        acc = na.omit(acc)
        tra = na.omit(tra)
        
        Test_set_stratified_accuracy <- mean(acc)
        Train_set_stratified_accuracy = mean(tra)
        trainsd <- sd(tra)
        testsd <- sd(acc)
        
        return (c(Test_set_stratified_accuracy,Train_set_stratified_accuracy, trainsd,testsd))
      }
      
      ###################### Calculate optimal split percentage ############################################
      # Calculate the optimal split percentage between .65 and .90 for predicting the final timepoint on the
      # exact previous one to be used in all splits for all timepoints. The criterion is the best accuracy  
      # as derived by the previous function with argument the different percentages.
      optimal_splitting_percentage <- function (infants){
        all_calculated_accuracies <- c()
        best_splitting_percentage <- c()
        for (percentage_to_calculate in 65:90){
          percentage_to_calculate = percentage_to_calculate/100
          telika <- multilogreg_stratified_One(infants = infants, optimal_splitting_percentage = percentage_to_calculate, timepoint = tail(timepoints_to_perform,1), end_timepoint = end_timepoint,meta_file = files,Ncomb = Ncomb,splitting_times = splitting_times)
          Test_set_stratified_accuracy = telika[1]
          Train_set_stratified_accuracy = telika[2]
          trainsd = telika[3]
          testsd = telika[4]
          all_calculated_accuracies= c(all_calculated_accuracies,Test_set_stratified_accuracy)
          best_splitting_percentage = c(best_splitting_percentage,percentage_to_calculate)
        }
        
        return (best_splitting_percentage[which.max(all_calculated_accuracies)])
      }
      
      # Assign the best percentage on a variable to use on the final calculations
      optimal_splitting_percentage <- optimal_splitting_percentage(infants = infants)
      
      ###################### Initializing variables Combinations ########################################
      Test_set_LOO_accuracy <- c()
      Train_set_LOO_accuracy <-c()
      Test_set_stratified_accuracy <- c()
      Train_set_stratified_accuracy <- c()
      trainsd <- c()
      testsd <- c()
      weighted_modeling_factors <- list()
      
      ###################### Calculate Accuracies on timepoints Combinations ############################
      
      # Loop to create a model from every timepoint to the final timepoint as selected from the first loop
      for (timepoint in timepoints_to_perform){
        
        LOOTeliko <- LOO_multilogreg_One(infants = infants, timepoint = timepoint,end_timepoint = end_timepoint,meta_file = files,Ncomb = Ncomb)
        
        Train_set_LOO_accuracy = c(Train_set_LOO_accuracy,LOOTeliko[1])
        Test_set_LOO_accuracy = c(Test_set_LOO_accuracy,LOOTeliko[2])
        
        MultiLogReg_results <- multilogreg_stratified_One(infants = infants, optimal_splitting_percentage = optimal_splitting_percentage,timepoint =  timepoint,meta_file =  files, end_timepoint = end_timepoint, Ncomb = Ncomb,splitting_times = splitting_times)
        
        
        Test_set_stratified_accuracy <- c(Test_set_stratified_accuracy,MultiLogReg_results[1])
        Train_set_stratified_accuracy <- c(Train_set_stratified_accuracy, MultiLogReg_results[2])
        trainsd <- c(trainsd, MultiLogReg_results[3])
        testsd <- c(testsd, MultiLogReg_results[4])
      }  
      
      
      ###################### Save accuracies on the lists ######################################
      
      # Select the effects that created the model
      effect = paste(effector_combinations[combination,],collapse = ' & ')
      
      # Adding to the lists in order to export the correct matrix
      Train_sets_LOO = rbind (Train_sets_LOO, c(effect,rep(0,(ncol(infants)-length(Train_set_LOO_accuracy)-1)),   Train_set_LOO_accuracy))
      Test_sets_LOO  = rbind (Test_sets_LOO,  c(effect,rep(0,(ncol(infants)-length(Test_set_LOO_accuracy)-1)),    Test_set_LOO_accuracy))
      Test_sets_stratified_split  = rbind (Test_sets_stratified_split,  c(effect,rep(0,(ncol(infants)-length(Test_set_stratified_accuracy)-1)), Test_set_stratified_accuracy))
      Train_sets_stratified_split = rbind (Train_sets_stratified_split, c(effect,rep(0,(ncol(infants)-length(Train_set_stratified_accuracy)-1)),Train_set_stratified_accuracy))
      
    }
  }
  
  
  
  
  
  ###################### Leave one out All #######################################
  # Calculate mean accuracy of prediction on final timepoint, via Leave One Out method of splitting Train-Test sets
  LOO_multilogreg <- function (infants,timepoint, end_timepoint){
    
    
    # Create the matrix with columns the state on the timepoint, the metadata and as a response the state on the last timepoint
    logreg <- as.data.frame(matrix(NA,nrow = nrow(infants),ncol = ncol(meta_file)))
    colnames(logreg)= c(paste('Cluster_at',end_timepoint,sep = '_'), timepoint , colnames(meta_file)[3:ncol(meta_file)])
    rownames(logreg)= rownames(infants)
    
    for (modeling_sample in (rownames(infants))){
      subTab = meta_file[meta_file[,1] %in% modeling_sample & meta_file[,2] == as.numeric(end_timepoint), 3:ncol(meta_file)]
      if (nrow(subTab) == 1){
        logreg[modeling_sample, colnames(meta_file)[3:ncol(meta_file)]] = subTab[1,]
      }
      logreg[modeling_sample, 2] = ifelse(test = length(infants[modeling_sample,timepoint]), yes = infants[modeling_sample,timepoint], no = NA)
      logreg[modeling_sample, 1] = ifelse(test = length(infants[modeling_sample,as.character(end_timepoint)]), yes = infants[modeling_sample,as.character(end_timepoint)], no = NA)
    }
    
    # Remove NAs
    logreg = logreg[complete.cases(logreg),]
    
    # Set the response variable
    fo = as.formula(paste(paste('Cluster_at',end_timepoint,sep = '_'),'~.',sep = ''))
    
    acc <- c()
    tra <- c()
    for (j in 1:nrow(logreg)){
      
      # Setting training and test sets by leaving one out
      trainset = logreg[1:nrow(logreg)!=j,]
      testset = logreg[j,]
      
      # Train the model
      prediction_model <- multinom(formula = fo , data = trainset ,censored = T, model = T)
      
      # Predict for trainset
      training_set_accuracies <- prediction_model %>% predict(trainset)
      # Predict for testset
      test_set_predictions <- prediction_model %>% predict(testset)
      # Calculate accuracies
      acc[j] = (sum(test_set_predictions == testset[,1])/ nrow(testset)) *100
      tra[j] = (sum(test_set_predictions == trainset[,1])/nrow(trainset))*100
      
    }
    
    return (c(mean(tra),mean(acc)))
  }
  
  ###################### Stratified Train/Test splits All ########################
  
  # Calculate accuracy of prediction on final timepoint, as a mean of 100 accuracies with different
  # Train-Test stratified splits. Stratfication is performed on all metadata categories 
  multilogreg_stratified <- function(infants,optimal_splitting_percentage,timepoint,times, end_timepoint,splitting_times){
    
    # Create the matrix with columns the state on the timepoint, the metadata and as a response the state on the last timepoint
    logreg <- as.data.frame(matrix(NA,nrow = nrow(infants),ncol = ncol(meta_file)))
    colnames(logreg)= c(paste('Cluster_at',end_timepoint,sep = '_'),timepoint,colnames(meta_file)[3:ncol(meta_file)])
    rownames(logreg)= rownames(infants)
    
    for (modeling_sample in (rownames(infants))){
      subTab = meta_file[meta_file[,1] %in% modeling_sample & meta_file[,2] == as.numeric(end_timepoint), 3:ncol(meta_file)]
      if (nrow(subTab) == 1){
        logreg[modeling_sample, colnames(meta_file)[3:ncol(meta_file)]] = subTab[1,]
      }
      logreg[modeling_sample, 2] = ifelse(test = length(infants[modeling_sample,timepoint]), yes = infants[modeling_sample,timepoint], no = NA)
      logreg[modeling_sample, 1] = ifelse(test = length(infants[modeling_sample,as.character(end_timepoint)]), yes = infants[modeling_sample,as.character(end_timepoint)], no = NA)
    }
    # Remove NAs
    logreg = logreg[complete.cases(logreg),]
    # Create different partitions to split train and test sets
    train_index <- createDataPartition(y = logreg[,paste('Cluster_at',end_timepoint,sep = '_')], p = optimal_splitting_percentage, list = F, times = splitting_times, groups = max(apply(X = logreg, MARGIN = 2, FUN = function(x){length(unique(x))})) )
    
    # Set the response variable
    fo = as.formula(paste(paste('Cluster_at',end_timepoint,sep = '_'),'~.',sep = ''))
    
    acc <- c()
    tra <- c()
    weighted_modeling_factors <- list()
    for (j in 1:ncol(train_index)){
      
      # Split train and test sets
      trainset = logreg[train_index[,j],]
      testset = logreg[-train_index[,j],]
      if (all(apply(X = logreg, MARGIN = 2, FUN = function(x){length(unique(x))})== apply(X = trainset, MARGIN = 2, FUN = function(x){length(unique(x))}))){
        
        # Create a model
        prediction_model <- multinom(formula = fo , data = trainset ,censored = T, model = T)
        # Predict for the trainset
        training_set_accuracies <- prediction_model %>% predict(trainset)
        # Predict for the testset
        test_set_predictions <- prediction_model %>% predict(testset)
        
        # Calculate accuracies
        acc[j] = (sum(test_set_predictions == testset[,1])/ nrow(testset)) *100
        tra[j] = (sum(training_set_accuracies == trainset[,1])/nrow(trainset))*100
      }
    }
    # Remove NAs
    acc = na.omit(acc)
    tra = na.omit(tra)
    
    Test_set_stratified_accuracy <- mean(acc)
    Train_set_stratified_accuracy = mean(tra)
    trainsd <- sd(tra)
    testsd <- sd(acc)
    
    return (c(Test_set_stratified_accuracy,Train_set_stratified_accuracy, trainsd,testsd))
  }
  
  ###################### Calculate optimal split percentage ######################
  # Calculate the optimal split percentage between .65 and .90 for predicting the final timepoint on the
  # exact previous one to be used in all splits for all timepoints. The criterion is the best accuracy  
  # as derived by the previous function with argument the different percentages.
  optimal_splitting_percentage <- function (infants){
    all_calculated_accuracies <- c()
    best_splitting_percentage <- c()
    for (percentage_to_calculate in 65:90){
      percentage_to_calculate = percentage_to_calculate/100
      telika <- multilogreg_stratified(end_timepoint = end_timepoint,infants = infants, optimal_splitting_percentage = percentage_to_calculate, timepoint = rev(colnames(infants))[2], splitting_times = splitting_times)
      Test_set_stratified_accuracy = telika[1]
      Train_set_stratified_accuracy = telika[2]
      trainsd = telika[3]
      testsd = telika[4]
      all_calculated_accuracies= c(all_calculated_accuracies,Test_set_stratified_accuracy)
      best_splitting_percentage = c(best_splitting_percentage,percentage_to_calculate)
    }
    
    return (best_splitting_percentage[which.max(all_calculated_accuracies)])
  }
  
  # Assign the best percentage on a variable to use on the final calculations
  optimal_splitting_percentage <- optimal_splitting_percentage(infants = infants)
  
  
  ###################### Initializing variables All #############################
  
  Test_set_LOO_accuracy <- c()
  Train_set_LOO_accuracy <-c()
  Test_set_stratified_accuracy <- c()
  Train_set_stratified_accuracy <- c()
  trainsd <- c()
  testsd <- c()
  weighted_modeling_factors <- list()
  
  
  ###################### Calculate Accuracies on timepoints All #################
  
  for (timepoint in timepoints_to_perform){
    
    
    LOOTeliko <- LOO_multilogreg(infants = infants, timepoint = timepoint,end_timepoint = end_timepoint)
    
    Train_set_LOO_accuracy = c(Train_set_LOO_accuracy,LOOTeliko[1])
    Test_set_LOO_accuracy = c(Test_set_LOO_accuracy,LOOTeliko[2])
    
    MultiLogReg_results <- multilogreg_stratified(infants, optimal_splitting_percentage, timepoint, times = 100,end_timepoint = end_timepoint,splitting_times = splitting_times)
    
    
    Test_set_stratified_accuracy <- c(Test_set_stratified_accuracy,MultiLogReg_results[1])
    Train_set_stratified_accuracy <- c(Train_set_stratified_accuracy, MultiLogReg_results[2])
    trainsd <- c(trainsd, MultiLogReg_results[3])
    testsd <- c(testsd, MultiLogReg_results[4])
  }  
  
  ###################### Save accuracies on the lists ######################################
  Train_sets_LOO = rbind (Train_sets_LOO,  c('All',rep(0,(ncol(infants)-length(Train_set_LOO_accuracy)-1)),Train_set_LOO_accuracy))
  Test_sets_LOO  = rbind (Test_sets_LOO,   c('All',rep(0,(ncol(infants)-length(Test_set_LOO_accuracy)-1)),Test_set_LOO_accuracy))
  Test_sets_stratified_split  = rbind (Test_sets_stratified_split,   c('All',rep(0,(ncol(infants)-length(Test_set_stratified_accuracy)-1)),Test_set_stratified_accuracy))
  Train_sets_stratified_split = rbind (Train_sets_stratified_split,  c('All',rep(0,(ncol(infants)-length(Train_set_stratified_accuracy)-1)),Train_set_stratified_accuracy))
  
  ###################### Barplots of Accuracies ############################################
  jpeg(filename = paste(output_dir,paste('Accuracies of Models on Timepoint', end_timepoint,sep = ' '),sep = '/'))
  barplot(rbind(Test_set_LOO_accuracy,Test_set_stratified_accuracy, Train_set_LOO_accuracy,Train_set_stratified_accuracy), beside = T, ylim = c(40,100), xpd = F ,names = timepoints_to_perform, col= c('paleturquoise4','paleturquoise3', 'darkolivegreen2', 'darkolivegreen3'))
  legend("topright", c("Test_sets_LOO","Test_sets_stratified_splitatified", 'Train_sets_LOO','Train_sets_stratified_splitatified'), fill=c('paleturquoise4','paleturquoise3', 'darkolivegreen2', 'darkolivegreen3'), cex = 1)
  abline(h = 40.3, col = "white", lwd = 2, lty = 2)
  
  dev.off()
  
  
}

###################### Calculate all possible combinations of metadata calculated ##########################
Npredictions = 0
for (r in 1:length(effectors)){
  Npredictions= Npredictions +(nrow(combinations(n = length(effectors),r = r,v = effectors,set = T,repeats.allowed = F)))
}
Npredictions = Npredictions +1

###################### Write files of All Calculated Accuracies ############################################
rownames(Train_sets_stratified_split) = Train_sets_stratified_split[,1]
rownames(Train_sets_LOO)= Train_sets_LOO[,1]
rownames(Test_sets_LOO) = Test_sets_LOO[,1]
rownames(Test_sets_stratified_split)= Test_sets_stratified_split[,1]
Train_sets_stratified_split = Train_sets_stratified_split[,2:ncol(Train_sets_stratified_split)]
Test_sets_stratified_split = Test_sets_stratified_split[,2:ncol(Test_sets_stratified_split)]
Train_sets_LOO =Train_sets_LOO[,2:ncol(Train_sets_LOO)]
Test_sets_LOO = Test_sets_LOO[,2:ncol(Test_sets_LOO)]

colnames(Train_sets_LOO) = c(paste('From timepoint ',rev(colnames(infants)[1:ncol(infants)-1]),sep = ''))
colnames(Test_sets_LOO)  = c(paste('From timepoint ',rev(colnames(infants)[1:ncol(infants)-1]),sep = ''))
colnames(Train_sets_stratified_split) = c(paste('From timepoint ',rev(colnames(infants)[1:ncol(infants)-1]),sep = ''))
colnames(Test_sets_stratified_split)  = c(paste('From timepoint ',rev(colnames(infants)[1:ncol(infants)-1]),sep = ''))
write.table(x = Train_sets_LOO, file = paste(output_dir,'All Accuracies of Training Sets LOO.csv', sep = '/') , col.names = F ,row.names = T)
write.table(x = Test_sets_LOO,  file = paste(output_dir,'All Accuracies of Test Sets LOO.csv',     sep = '/') , col.names = F,row.names = T)
write.table(x = Train_sets_stratified_split, file = paste(output_dir,'All Accuracies of Training Sets Splitted.csv', sep = '/') , col.names = F, row.names = T)
write.table(x = Test_sets_stratified_split,  file = paste(output_dir,'All Accuracies of Test Sets Splitted.csv', sep = '/')  , col.names = F,row.names = T)


###################### Calculate random estimators performance ##########################################
random_estimator = 100/rev(apply (X = infants,MARGIN = 2,FUN = function(x){max(x,na.rm = T)}))
###################### Write files of Best Accuracies Test_sets_LOO ###########################################
Test_sets_stratified_split

TotimepointIndeces = c(0:ncol(Test_sets_LOO))*Npredictions
maxaccuracies <- c()
metadata <- c()
for (i in 1:(length(TotimepointIndeces)-1)){
  
  maxaccuracies[i]   = max(Test_sets_LOO[TotimepointIndeces[i]:TotimepointIndeces[i+1],i])
  metadata[i] = rownames(Test_sets_LOO)[which.max(Test_sets_LOO[TotimepointIndeces[i]:TotimepointIndeces[i+1],i])]
  
}
write.csv(x = rbind(paste(paste('Maximum Accuracy for Timepoint', rev(colnames(infants)[2:ncol(infants)])), c(rev(colnames(infants))[2:ncol(infants)]), sep = ' from Timepoint '),maxaccuracies,metadata,random_estimator[1:length(random_estimator)-1]),file = paste(output_dir,'Maximum_Accuracies_of_Test_sets_LOO.csv',sep = '/'),row.names = F)

###################### Write files of Best Accuracies TestSplits ########################################

maxaccuracies <- c()
metadata <- c()
for (i in 1:(length(TotimepointIndeces)-1)){
  
  maxaccuracies[i]   = max(Test_sets_stratified_split[TotimepointIndeces[i]:TotimepointIndeces[i+1],i])
  metadata[i] = rownames(Test_sets_stratified_split)[which.max(Test_sets_stratified_split[TotimepointIndeces[i]:TotimepointIndeces[i+1],i])]
  
}

write.csv(x = rbind(paste(paste('Maximum Accuracy for Timepoint', rev(colnames(infants)[2:ncol(infants)])), c(rev(colnames(infants))[2:ncol(infants)]), sep = ' from Timepoint '),maxaccuracies,metadata,random_estimator[1:length(random_estimator)-1]),file = paste(output_dir,'Maximum_Accuracies_of_Stratified_Tests.csv',sep = '/'),row.names = F)

###################### Write files of Best Accuracies Train_sets_LOO ########################################

maxaccuracies <- c()
metadata <- c()
for (i in 1:(length(TotimepointIndeces)-1)){
  
  maxaccuracies[i]   = max(Train_sets_LOO[TotimepointIndeces[i]:TotimepointIndeces[i+1],i])
  metadata[i] = rownames(Train_sets_LOO)[which.max(Train_sets_LOO[TotimepointIndeces[i]:TotimepointIndeces[i+1],i])]
  
}

write.csv(x = rbind(paste(paste('Maximum Accuracy for Timepoint', rev(colnames(infants)[2:ncol(infants)])), c(rev(colnames(infants))[2:ncol(infants)]), sep = ' from Timepoint '),maxaccuracies,metadata,random_estimator[1:length(random_estimator)-1]),file = paste(output_dir,'Maximum_Accuracies_of_Train_sets_LOO.csv',sep = '/'),row.names = F)

###################### Write files of Best Accuracies TrainSplits ########################################

maxaccuracies <- c()
metadata <- c()

for (i in 1:(length(TotimepointIndeces)-1)){
  
  maxaccuracies[i]   = max(Train_sets_stratified_split[TotimepointIndeces[i]:TotimepointIndeces[i+1],i])
  metadata[i] = rownames(Train_sets_stratified_split)[which.max(Train_sets_stratified_split[TotimepointIndeces[i]:TotimepointIndeces[i+1],i])]
  
}

write.csv(x = rbind(paste(paste('Maximum Accuracy for Timepoint', rev(colnames(infants)[2:ncol(infants)])), c(rev(colnames(infants))[2:ncol(infants)]), sep = ' from Timepoint '),maxaccuracies,metadata,random_estimator[1:length(random_estimator)-1]),file = paste(output_dir,'Maximum_Accuracies_of_TrainSplits.csv',sep = '/'),row.names = F)

###################### Perform Chi-square analysis to check metadata influence on T0 #####################

# Create a matrix to check for the effect of the metadata on the first timepoint (T0)
chimatrix = matrix (NA, ncol = length(effectors)+1, nrow = nrow(infants))

# Assign the rownames to the matrix
rownames(chimatrix) = rownames(infants)

# Assign the values on the matrix
chimatrix[,1] = infants[,1]

for (chirow in rownames(infants)){
  
  subTab = as.matrix(meta_file[meta_file[,'Sample']==chirow & meta_file[,'Timepoint']==1,3:ncol(meta_file)])
  
  if (length(subTab)){
    chimatrix[chirow,2:ncol(chimatrix)] = subTab
  }
}
chimatrix = chimatrix[complete.cases(chimatrix),]

Nuniq = length(unique(chimatrix[,1]))

Significant_metadata = c()
for (metadatum in 2:ncol(chimatrix)){
  
  chi_square_value = 0
  dof = 0
  
  for (secondstate in unique(chimatrix[,metadatum])){
    
    expected = sum(chimatrix[,metadatum]==secondstate) * (1/length(unique(chimatrix[,1])))
    
    for (firstcluster in unique(chimatrix[,1])){
      dof = dof + 1
      observed = sum(chimatrix[chimatrix[,metadatum]== secondstate,1]==firstcluster)
      chi_square_value = chi_square_value + ((observed- expected)**2)/expected 
    }
  }
  
  Null_hypothesis = chi_square_value <= qchisq(.95, df = (dof-1),lower.tail = T)
  Significant_metadata[metadatum] = Null_hypothesis
}

Significant_metadata = Significant_metadata[2:length(Significant_metadata)]

if (any(Significant_metadata)){
  Significant_metadata = colnames(meta_file)[3:ncol(meta_file)][Significant_metadata]
}




print (' Analysis Completed ')
print ('_____________________________________________')
print ('Results are saved in the preselected folders ')
print ('_____________________________________________')



} else {
  print (paste('Cronos has already run with the exact same parameters. The output files are stored in',sub(pattern = './',replacement = '',x = directory),sep = ' '))
}
