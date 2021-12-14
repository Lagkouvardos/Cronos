##########################################################################################
############ ON THIS SECTION YOU CAN SELECT YOUR DESIRED PARAMETERS ######################
##########################################################################################

########### PLEASE FOLLOW THE INSTRUCTIONS CAREFULLY #####################################

# Please set the directory of the script as the working folder (e.g D:/studyname/Data/Chronos/)
# Note: the path is denoted by forward slash "/".
working_directory = "~/Working_Chronos/Cronos_Final/" #<--- CHANGE ACCORDINGLY !!!   "~/Working_Chronos/Cronos_Final/"

# Please give the file name of the normalized OTU-table without taxonomic classification
input_otu = "SOTUs-Table.tab"           #<--- CHANGE ACCORDINGLY !!!
# Please give the name of the meta-file that contains individual sample information
input_meta = "Mappings_of_ot&24&St.csv"                #<--- CHANGE ACCORDINGLY !!!
# Please give the name of the phylogenetic tree constructed from the OTU sequences
input_tree = "SOTUs-NJTree-All.tre"         #<--- CHANGE ACCORDINGLY !!!


# Please specify if the file contains external data. One example is when analyzing infant data
# but have an OTU table for adults as well for a reference. 
# If it contains external timepoint, please specify the name of the column 
# (i.e. timepoint name) where they are saved
# If it doesn't, leave it blank (like this: '' ).
External_Reference_Point = 'MM'              #<--- CHANGE ACCORDINGLY


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
# Greater numbers lead to more time needed (default is 10)
# Must be over 1
splitting_times = 2                           # <---- CHANGE ACCORDINGLY

# Please select the action Cronos should do if the same parameters are already used for analysis before
# It can be either 'Continue' or 'Stop'.
action = 'Continue'                           # <---- CHANGE ACCORDINGLY


########### Now that you selected the parameters you can select everything (Ctrl+A) and press run (Ctrl+Enter)  ##########

##########################################################################################################################
###################### PLEASE DO NOT CHANGE ANYTHING BEYOND THIS POINT ###################################################
##########################################################################################################################

# Here we set the working directory as you selected in order to find your files and R scripts
# Please do not interfere with the next line

setwd(working_directory)

##########################################################################################################################
#################################### CREATING OUTPUT DIRECTORY ###########################################################
##########################################################################################################################

# Setting the name of the directory where Cronos outputs are stored 
date_of_run = unlist(strsplit(as.character(date()), split = ' '))
date_of_run = date_of_run[nchar(date_of_run)>0]  
date_of_run = paste(date_of_run,collapse = '_')
date_of_run= gsub(pattern = ':', replacement = '_', x = date_of_run)
output_dir = paste('Cronos',date_of_run,sep ='_')
# Create a log file with all the parameters used
parameters = matrix('',ncol = 1 ,nrow = 6)
rownames(parameters) = c('input_meta','input_tree','input_otu','taxonomic_level','External_Reference_Point','splitting_times')
parameters[,1]= c(input_meta,input_tree,input_otu,taxonomic_level,External_Reference_Point,splitting_times)

# Checking whether Cronos was run before with the same parameters
new_run = T
for (directory in list.dirs()[2:length(list.dirs())]){
  previous_runs = list.files(path = directory)
  existing_runs = read.csv(paste(directory,previous_runs[grepl(pattern = '_log.tab', x = previous_runs , ignore.case = F)], sep = '/'))
  if (all(existing_runs[,2] == parameters[,1])){
    new_run = F
  }
}

if (new_run==T || (new_run == F & action =='Continue')){
  
  # Create the directory where Cronos outputs will be stored
  dir.create(paste(working_directory,output_dir,sep = '/'),showWarnings = F)
  # Write file with the parameters on this run
  
  write.table(x = parameters,file = paste(output_dir,'Cronos_log.tab',sep = '/'), sep = "\t",col.names =F, row.names = TRUE,quote = FALSE)
  

  ##################################### SECTION ############################################################################
  #################################### CLUSTERING ##########################################################################
  ##########################################################################################################################
  
  ############################ Reading the files ###########################################
  meta_file <- read.table (file = input_meta, check.names = FALSE, header = TRUE, dec = ".", sep = "\t", row.names = 1, comment.char = "", stringsAsFactors = F)
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
  
  packages <-c("ade4","dplyr","GUniFrac","phangorn","cluster","fpc","markovchain", 'spgs','caret','nnet','gtools', 'mclust','igraph', 'network','ggplot2','reshape2','easyalluvial','ggrepel', 'vegan') 
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
  
  # Vectors to save medoid-relative information
  medoids_plot <- c()
  medoid_names <- c()
  number_samples <- c()
  time <- c()
  clusters_medoids <- c()
  
  # Calculate the UniFrac distance matrix for comparing microbial communities
  for (name in names(timepoint_list)){
    
    # Create temporary files
    temp_clusters_medoids <- c()
    temp_medoid_names <- c()
    temp_number_samples <- c()
    temp_time <- c()
    
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
        calinski_harabasz_values= c(calinski_harabasz_values,(cluster.stats(unifract_dist,clusteringP)[["ch"]]))
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
    
    dir.create(paste(output_dir,"Calinski-Harabasz",sep = '/'),showWarnings = F)
    
    # Export plot with the Calinski-Harabasz scores for each k
    jpeg(filename =paste(paste(output_dir,"Calinski-Harabasz",sep = '/'), paste(paste('Calinski-Harabasz_index',name,sep = ' of Timepoint '),'jpeg',sep = ' .'),sep = '/'),width = 800,height=842)
    plot(2:9,calinski_harabasz_values, type="h", xlab="k clusters", ylab="CH index",main='Calinski-Harabasz scores of different #Clusters')
    legend("topright", c("PAM"), fill=c("black"))
    dev.off()
    
    # Perform clustering
    clustering_results <- PAM_clustering(unifract_dist = unifract_dist ,k =best_k)
    clusters = clustering_results[[1]]
    medoids_temp = rep(NA,9-length(clustering_results[[2]]))
    medoids[,name] = c(clustering_results[[2]], medoids_temp)
    
    # Store the information about the medoids, the timepoints and the clusters
    if (name!="ot"){
      if (name!=External_Reference_Point){
        temp_medoids_plot = c(clustering_results[[2]])
        for (i in 1:best_k){
          temp_medoid_names <- c(temp_medoid_names, paste0("TP:",name," Cl:",i))
          temp_number_samples <- c(temp_number_samples,length(which(clusters==i))/length(clusters))
          temp_time <- c(temp_time,name)
          temp_clusters_medoids <- c(temp_clusters_medoids,paste("Cluster",i))
        }
      } else {
        temp_medoids_plot = c( clustering_results[[2]])
        for (i in 1:best_k){
          temp_medoid_names <- c( temp_medoid_names,paste0(External_Reference_Point," Cl:",i))
          temp_number_samples <- c(temp_number_samples,length(which(clusters==i))/length(clusters))
          temp_time <- c(temp_time,name)
          temp_clusters_medoids <- c(temp_clusters_medoids,paste("Cluster",i))
          
        }
      }
    }
    
    #  medoids = c(medoids, clustering_results[[2]])
    avg_width = clustering_results[[3]]
    # Assing the samples to the estimated clusters
    samples_on_clusters[meta_file[row.names(unifract_dist),'Sample'],name]<- clusters
    samples_on_clusters[samples_on_clusters==0] <- NA
    
    
    # Function for the MDS Plots
    sclass = function(distance,groups,colors) {
      
      # Performing classical multidimensional scaling
      mds<- cmdscale(distance,eig=T, x.ret=T)
      # Computing the variance of each axis
      mds.variation<- round(mds$eig/sum(mds$eig)*100,1) 
      all_groups_comp = groups[!is.na(groups)]
      
      # Performing analysis of variance using the distance matrix(only if individuals=NULL)
      if(nlevels(groups) > 1) {
        # PERMANOVA test
        adonis = adonis(distance ~ all_groups_comp)
        # Create Subtitle
        sub = paste("MDS plot of Microbial Profiles\n(p-value ",adonis[[1]][6][[1]][1],")",sep="")
      } else {
        sub=c("")
      }
      # Find the factors of all_groups_comp vector
      all_groups_comp = factor(all_groups_comp,levels(all_groups_comp)[unique(all_groups_comp)])
      
      # Display the MDS plot (Multidimensional Scaling plot)
      s.class(
        mds$points, col = colors, cpoint =
          2, fac = all_groups_comp, grid = T, sub = sub)
      graphics:: title (main=paste0("MDS graph " ), xlab = paste("MDS1:",mds.variation[1],"%"), ylab=paste("MDS2:",mds.variation[2],"%"))
    }
    
    
    dir.create(paste(output_dir,"MDS plots",sep = '/'),showWarnings = F)
    
    # Export MDS plot of the samples on the timepoint
    pdf(paste(paste(output_dir,"MDS plots",sep = '/'),paste( name,'MDS Plot of Timepoint.pdf',sep=' '),sep = '/'))
    sclass(unifract_dist,as.factor(clusters), colours_ploting[1:nlevels(as.factor(clusters))])
    dev.off()
    
    jpeg(filename = paste(paste(output_dir,"MDS plots",sep = '/'),paste( name,'MDS Plot of Timepoint.jpeg',sep=' '),sep = '/'),width = 800,height=842)
    sclass(unifract_dist,as.factor(clusters), colours_ploting[1:nlevels(as.factor(clusters))])
    dev.off()
    
    #############Gaussian Mixture Model Based Test to evaluate whether the dataset forms groups ###############################################
    
    gmm_testing = Mclust(data = as.dist(unifract_dist),G = c(1,max(samples_on_clusters[,name], na.rm = T)),  modelNames = c("EII","VII","EEI","EVI","VEI","VVI") , verbose = F)
    
    BIC_list[[name]] = gmm_testing$BIC
    if (gmm_testing$G == 1){
      samples_on_clusters[,name] = rep(1,nrow(samples_on_clusters))
      
      # Store the information about the medoids, the timepoints and the clusters
      if (name!="ot"){
        if (name!=External_Reference_Point){
          clustering_results = PAM_clustering(unifract_dist = unifract_dist, k = 1)
          temp_medoids_plot = c( clustering_results[[2]])
          
          temp_medoid_names <- c( paste0("TP:",name," Cl:",1))
          temp_number_samples <- c(length(which(clusters==1))/length(clusters))
          temp_time <- name
          temp_clusters_medoids <- c(paste("Cluster",1))
          
        } else {
          clustering_results = PAM_clustering(unifract_dist = unifract_dist, k = 1)
          temp_medoids_plot = c(clustering_results[[2]])
          
          medoid_names <- c(medoid_names, paste0(External_Reference_Point," Cl:",1))
          temp_number_samples <- c(length(which(clusters==1))/length(clusters))
          temp_time <- name
          temp_clusters_medoids <- c(paste("Cluster",1))
          
          
        }
        
      }
    }
    
    # Store the information about the medoids, the timepoints and the clusters
    if (name!="ot"){
      medoids_plot = c(medoids_plot, temp_medoids_plot)
      medoid_names <- c(medoid_names, temp_medoid_names)
      number_samples <- c(number_samples,temp_number_samples)
      time <- c(time,temp_time)
      clusters_medoids <- c(clusters_medoids,temp_clusters_medoids)
    }
    
  }
  
  
  ######################################### SECTION ########################################################################
  #################################### TRANSITION ANALYSIS #################################################################
  ##########################################################################################################################
  
  ############ Markovian property check ################################################
  
  state_calculator = function(dataset_full,t0,t2){
    states <- c()
    for (j in 1:nrow(unique(dataset_full[,t0:t2]))){
      if (!any(is.na(unique(dataset_full[,t0:t2])[j,]))){
        states[[j]] <- as.vector (unique(dataset_full[,t0:t2])[j,])
      }
    }
    return (states)  
  }  
  
  counting_states = function (dataset_full,t0,t2,states) {
    counter <- c(rep(0,length(states)))
    for (i in 1:nrow(dataset_full)){
      if (!any(is.na(unique(dataset_full[i,t0:t2])))){
        for (j in 1:length(states)){
          if (all(dataset_full[i,t0:t2]==states[[j]])){
            counter[j]= counter[j]+1
          } 
        }
      }
    }
    return (counter)
  }
  
  if (nchar(External_Reference_Point)>0){
    tps = as.character(sort(as.numeric(colnames(samples_on_clusters)[colnames(samples_on_clusters)!='ot' & colnames(samples_on_clusters)!=External_Reference_Point]),decreasing = F))
  } else {
    tps = as.character(sort(as.numeric(colnames(samples_on_clusters[colnames(samples_on_clusters)!='ot'])),decreasing = F))
  }
  
  
  dataset_full = samples_on_clusters[,tps]
  dataset_full = dataset_full[!apply(dataset_full,1,function(x){all(is.na(x))}),]
  dataset_full_on_clusters <- dataset_full
  independed_clusters <- cumsum(as.vector(apply(X = dataset_full,MARGIN = 2,FUN = function(x){max(x,na.rm = T)})))
  
  
  ######## Auxiliary Matrix necessary for the  Transition Plot ##########
  
  # Calculate the number of samples per cluster. This is necessary for the plots
  samples_per_cluster <- data.frame()
  for (i in 1:ncol(dataset_full_on_clusters)){
    for (j in 1:max(dataset_full_on_clusters[,i],na.rm = T)){
      samples_per_cluster[j,i] <- length(which(c(dataset_full_on_clusters[,i])==j))
    }
  }
  
  rownames(samples_per_cluster) <- paste("Cluster",1:nrow(samples_per_cluster))
  
  
  # Replace NAs with zero
  samples_per_cluster[is.na(samples_per_cluster)] <- 0
  
  # Convert absolute values to ratios
  samples_per_cluster <-  t(t(samples_per_cluster) / nrow(dataset_full_on_clusters))
  
  # Unpivot samples_per_cluster from wide to long format using the clusters as the identifier
  samples_per_cluster <- data.frame(rownames(samples_per_cluster),samples_per_cluster)
  colnames(samples_per_cluster) <- c("Clusters",colnames(dataset_full_on_clusters))
  samples_per_cluster <-  melt(samples_per_cluster, cols = 2:ncol(samples_per_cluster),id.vars ="Clusters")
  colnames(samples_per_cluster) <- c("Clusters","Timepoints","Samples")
  
  # Convert Timepoints and Cluster columns into factors
  samples_per_cluster$Timepoints <- factor(samples_per_cluster$Timepoints,levels=colnames(dataset_full_on_clusters))
  samples_per_cluster$Clusters <- factor(samples_per_cluster$Clusters,levels=unique(samples_per_cluster$Clusters))
  
  # Colours of the bubble plot
  bubble_color <- colours_ploting[1:nlevels(as.factor(samples_per_cluster[,"Timepoints"]))]
  
  
  
  ######## Auxiliary Matrix necessary for the  Medoids MDS Plot ##########
  
  # Dataframe with the medoids
  medoids_plot <- data.frame(medoid_names,medoids_plot,number_samples)
  
  # Calculate the distances between the medoids
  unifracs_plot <- GUniFrac( otu_file[medoids_plot[,2],] ,tree = rooted_tree, alpha = c(0.0,0.5,1.0))$unifracs
  
  # Weight on abundant lineages so  the distance is not dominated by highly abundant lineages with 0.5 having the best power
  unifract_dist_plot <- unifracs_plot[, , "d_0.5"]
  
  # Perform MDS scaling
  mds<- cmdscale(unifract_dist_plot,eig=T, x.ret=T)
  # Calculate the axes variance
  mds.variation<- round(mds$eig/sum(mds$eig)*100,1) 
  # Coordinates of the points
  mds.points<- mds$points 
  
  # Dataframe with the coordinates and information about the medoids
  mdsdata<- data.frame(Sample=rownames(mds.points),X=mds.points[,1] ,Y=mds.points[,2],Samples=number_samples,Points=medoid_names,Time=time, Clusters=clusters_medoids)
  
  
  
  ######## Auxiliary Matrix necessary for the  Alluvial Plot ##########
  
  # Prepare the data for the alluvial plot
  alluvial_matrix <- dataset_full_on_clusters
  
  # Replace NAs with "NA"
  alluvial_matrix[is.na(alluvial_matrix)] <- "NA"
  
  # Unpivot alluvial_matrix from wide to long format using the IDs as the identifier
  alluvial_matrix <- data.frame(rownames(alluvial_matrix),alluvial_matrix)
  colnames(alluvial_matrix) <- c("IDs",colnames(dataset_full_on_clusters))
  alluvial_plot <- melt(alluvial_matrix, id.vars ="IDs")
  colnames(alluvial_plot) <- c("Samples","Timepoints","Clusters")
  
  # Convert Timepoints column into factors
  alluvial_plot[,"Timepoints"] <- factor( alluvial_plot[,"Timepoints"], levels=colnames(alluvial_matrix)[2:ncol(alluvial_matrix)])
  # Convert Timepoints column into characters
  alluvial_plot[,"Clusters"] <- as.character(alluvial_plot[,"Clusters"])
  
  # Colours of the Alluvial plot
  alluvial_color <- colours_ploting[1:nlevels(as.factor(alluvial_plot[,"Clusters"]))]
  
  
  
  
  
  for (i in 2:ncol(dataset_full_on_clusters)){
    for (j in 1:max(dataset_full_on_clusters[,i],na.rm = T)){
      dataset_full_on_clusters[dataset_full_on_clusters[,i]==j,i] <- independed_clusters[i-1]+j
    }
  }
  
  Markov_Test_X2 <- function(dataset_full){
    Ss <- c()
    dof <- c()
    for (t0 in 1:(ncol(dataset_full)-2)){
      t1 <- t0+1
      t2 <- t0+2
      ijk_states <- state_calculator(dataset_full = dataset_full, t0 = t0 , t2 = t2 )
      ijk_counter <- counting_states(dataset_full = dataset_full ,t0 = t0, t2 = t2, states = ijk_states)
      
      jk_states <-lapply(ijk_states, function(x){x[2:3]})
      jk_counter <- counting_states(dataset_full = dataset_full ,t0 = t1,t2 = t2, states = jk_states)
      ij_states <- lapply(ijk_states, function(x){x[1:2]})
      ij_counter <- counting_states(dataset_full = dataset_full ,t0 = t0, t2 = t1,states = ij_states)
      
      p_jk <- jk_counter/sum(jk_counter)
      S =0
      for (i in 1:length(jk_counter)){
        S = S + (ijk_counter[i]-ij_counter[i]*p_jk[[i]])**2/ij_counter[i]*p_jk[[i]]
      }
      degrees_of_freedom = length(unique(dataset_full[,t0])) - length(unique(ij_states)) + length(unique(ijk_states)) +1  
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
  
  Markovian_Check <- Markov_Test_X2(dataset_full = dataset_full)
  Markovian_property= all(Markovian_Check[1,]==1)
  
  ############ Claim transition matrix ################################################
  
  # Vectors with the initial and final positions of each timepoint
  x_initial <- c()
  x_final <- c()
  y_initial <- c()
  y_final <- c()
  transitions <- c()
  
  matrix_of_transitions <- matrix(0, nrow = max(dataset_full_on_clusters, na.rm = T), ncol = max(dataset_full_on_clusters,na.rm = T))
  for (i in 1:(ncol(dataset_full_on_clusters)-1)){
    for (j in min(dataset_full_on_clusters[,i],na.rm = T):max(dataset_full_on_clusters[,i],na.rm = T)){
      for (k in min(dataset_full_on_clusters[,i+1],na.rm = T):max(dataset_full_on_clusters[,i+1],na.rm = T)){
        
        matrix_of_transitions[j,k] = sum(dataset_full_on_clusters[dataset_full_on_clusters[,i]==j,(i+1)]==k,na.rm = T)/sum(dataset_full_on_clusters[,i]==j,na.rm = T)
        if (matrix_of_transitions[j,k]>0){
          x_initial <- c(x_initial,colnames(dataset_full_on_clusters)[i])
          y_initial <- c(y_initial,paste("Cluster",(j-(min(dataset_full_on_clusters[,i],na.rm = T)-1))))
          x_final <- c(x_final,colnames(dataset_full_on_clusters)[i+1])
          y_final <- c(y_final,paste("Cluster",(k-(min(dataset_full_on_clusters[,i+1],na.rm = T)-1))))
          transitions <- c(transitions,matrix_of_transitions[j,k])
        }
        
      }
    }
  }
  
  # Vectors with the initial and final positions of each timepoint (Transition matrix for the MDS0transition plot)
  x_initial_mds <- c()
  x_final_mds <- c()
  y_initial_mds <- c()
  y_final_mds <- c()
  
  for(i in 1:length(x_initial)) {
    row <- which(mdsdata[,"Time"]==x_initial[i] & mdsdata[,"Clusters"]==y_initial[i])
    row2 <- which(mdsdata[,"Time"]==x_final[i] & mdsdata[,"Clusters"]==y_final[i])
    
    x_initial_mds <- c(x_initial_mds,mdsdata[row,2])
    y_initial_mds <- c(y_initial_mds,mdsdata[row,3])
    x_final_mds <- c(x_final_mds,mdsdata[row2,2])
    y_final_mds <- c(y_final_mds,mdsdata[row2,3])
  }
  
  # An alternative transition matrix. This is necessary for the MDS-transition plot
  mds_transition_matrix <- data.frame(x1 = x_initial_mds, x2=x_final_mds, y1 = y_initial_mds , y2 = y_final_mds)
  
  # An alternative transition matrix. This is necessary for the bubble plot
  bubble_transition_matrix <- data.frame(x1 = x_initial, x2=x_final, y1 = y_initial , y2 = y_final)
  
  # Convert Timepoints into factors
  if (nchar(External_Reference_Point)>0) {  
    mdsdata$Time <- factor(mdsdata$Time,levels=c(colnames(dataset_full_on_clusters),External_Reference_Point))
  } else {
    mdsdata$Time <- factor(mdsdata$Time,levels=colnames(dataset_full_on_clusters))
  }
  
  # Colours of the MDS transition plot
  mds_color <- colours_ploting[1:nlevels(as.factor(mdsdata$Time))]
  
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
  samples_on_clusters = samples_on_clusters[,c(as.character(sort(as.numeric(colnames(samples_on_clusters)[colnames(samples_on_clusters)!=External_Reference_Point]),decreasing = F,na.last = T)),External_Reference_Point)]
  
  
  dataset_full<- samples_on_clusters[is.na(samples_on_clusters[,External_Reference_Point]),1:ncol(samples_on_clusters)-1]
  
  ############ Write comma delimited files with outputs   ####################################
  
  colnames(taxa_clusters)<-lapply(X = colnames(taxa_clusters), FUN = function(x){paste('Timepoint',x,sep = ' ')})
  
  write.table(x = samples_on_clusters, file = paste(output_dir, "Samples_in_Timepoint-specific_Clusters.tab",sep = '/'), sep = "\t",col.names =NA, row.names = TRUE,quote = FALSE)
  write.table(x = taxa_clusters,       file = paste(output_dir, "Taxonomic_profile_of_clusters.tab", sep = '/'), sep = "\t",col.names =NA, row.names = TRUE,quote = FALSE)
  

  No_clusters_per_timepoint = apply(dataset_full, 2,function(x){max(x,na.rm = T)})
  transition_names = c()
  for (j in 1:length(No_clusters_per_timepoint)){
    for (i in 1:(No_clusters_per_timepoint[j])){
      transition_names= c(transition_names,paste( paste('Timepoint',colnames(dataset_full)[j],sep = ' '), paste( 'Cluster',i,sep = ' '),sep = ' '))
    }
  }
  rownames(matrix_of_transitions)= paste('From ', transition_names,sep = '')
  
  write.table(x = matrix_of_transitions,file = paste(output_dir,'Transition_Matrix.tab',sep = '/'), sep = "\t",col.names =NA, row.names = TRUE,quote = FALSE)
  
  
  
  ############ Plotting the transitions ####################################
  
  # Create a new directory where the transition plots will be placed
  dir.create(paste(output_dir,"Transitions plots",sep = '/'),showWarnings = F)
  
  # Transition-Bubble Plot
  bubble_plot <- ggplot(samples_per_cluster, aes(y = Clusters, x = Timepoints)) + 
    geom_point(aes(size = Samples, fill = Timepoints), alpha = 0.75, shape =21 ,show.legend = F) + 
    scale_y_reverse()+
    scale_size_continuous(limits = c(0.000001, 1), range = c(1,60), breaks = c(seq(0.000001,1+0.000001,by=0.01))) +
    labs( x= "", y = "" ,fill = "")  + 
    scale_y_discrete(limits=rev)+
    scale_x_discrete(position = "top") +
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2,color=transitions), data = bubble_transition_matrix,size=transitions*5)+
    scale_color_gradient(low = "grey90", high = "black", n.breaks=5)+
    geom_point(aes(size = Samples, fill = Timepoints), alpha = 1, shape =21 ,show.legend = F) + 
    theme( 
      legend.key=element_blank(),
      axis.text.x = element_text(colour = "black", size = 12, face = "bold"),
      axis.text.y = element_text(colour = "black", face = "bold", size = 11),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.background = element_blank(),
      legend.position = "right"
      ,plot.title = element_text(hjust = 0.5,size=20, face = "bold")) +
    ggtitle("Transition Plot")+
    scale_fill_manual(values = bubble_color, guide = "none") 
  
  
  # Print the transition-Bubble Plot 
  suppressWarnings(ggsave(filename = paste(paste(output_dir,"Transitions plots",sep = '/'),'Transition Plot.pdf',sep = '/'),bubble_plot,width = 400,height=300,unit="mm"))
  jpeg(filename = paste(paste(output_dir,"Transitions plots",sep = '/'),'Transition Plot.jpeg',sep = '/'),width = 1200 ,height=842)
  suppressWarnings(print(bubble_plot))
  dev.off()
  
  
  # Alluvial-Flow Chart
  alluvial <- alluvial_long(alluvial_plot , key = Timepoints , value = Clusters , id = Samples , fill_by = 'value' , col_vector_flow =  alluvial_color , col_vector_value = alluvial_color)+
    xlab("Timepoints")+
    ylab("Clusters")+
    ggtitle("Alluvial Plot")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5,size=20, face = "bold"),axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  # Print the alluvial plot
  ggsave(filename = paste(paste(output_dir,"Transitions plots",sep = '/'),'Alluvial Plot.pdf',sep = '/'),alluvial,width = 400,height=300,unit="mm")
  jpeg(filename = paste(paste(output_dir,"Transitions plots",sep = '/'),'Alluvial Plot.jpeg',sep = '/'),width = 842 ,height=842)
  show(alluvial)
  dev.off()
  
  
  # Medoids MDS Plot
  medoids_mds <-  ggplot(mdsdata,aes(x=X,y=Y))+
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2,color=transitions), data = mds_transition_matrix,arrow = arrow(length = unit(0.02, "npc")),size=transitions*2)+
    scale_color_gradient(low = "grey90", high = "black", n.breaks=5)+
    scale_size_continuous(limits = c(0.000001, 1), range = c(1,20), breaks = c(seq(0.000001,1+0.000001,by=0.01))) +
    geom_point(aes(size = Samples,fill=Time),shape=21,alpha=0.4,show.legend = F)+
    geom_text_repel(aes(label=Points),size=4)+
    xlab(paste("MDS1:",mds.variation[1],"%"))+
    ylab(paste("MDS2:",mds.variation[2],"%"))+
    ggtitle("Timepoints Medoids")+
    coord_fixed()+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_fill_manual(values = mds_color, guide = "none") 
  
  
  # Print the medoids_mds
  ggsave(filename = paste(paste(output_dir,"Transitions plots",sep = '/'),'Medoids MDS Plot.pdf',sep = '/'),medoids_mds,width = 500,height=300,unit="mm")
  jpeg(filename = paste(paste(output_dir,"Transitions plots",sep = '/'),'Medoids MDS Plot.jpeg',sep = '/'),width = 1200 ,height=842)
  show(medoids_mds)
  dev.off()
  
  ##################################### SECTION ############################################################################
  #################################### MODELING ############################################################################
  ##########################################################################################################################
  
  heatmap_matrix_LOO <- data.frame(matrix(NA,nrow = ncol(dataset_full_on_clusters), ncol=ncol(dataset_full_on_clusters)))
  diag(heatmap_matrix_LOO) <- 100
  
  heatmap_matrix_statified <- data.frame(matrix(NA,nrow = ncol(dataset_full_on_clusters) , ncol=ncol(dataset_full_on_clusters)))
  diag(heatmap_matrix_statified) <- 100
  
  ###################### Initializing Lists to save metrics of modeling ############################
  # Initializing lists to write accuracies of models
  Train_sets_LOO <- c()
  Test_sets_LOO  <- c()
  Train_sets_stratified_split <- c()
  Test_sets_stratified_split  <- c()
  Train_LOO_overtime_accuracy <- c()
  Test_LOO_overtime_accuracy <- c()
  Train_stratified_overtime_accuracy <- c()
  Test_stratified_overtime_accuracy <- c()
  # Naming the effectors to be imputed on the model
  effectors = colnames(meta_file)[3:ncol(meta_file)]
  
  ###################### Perform Tasks on all Timepoints ###############################
  
  # Loop to create a model for each timepoint on the dataset
  for (end_timepoint in head(as.numeric(rev(colnames(dataset_full))),n = (ncol(dataset_full))-1)){
    
    
    # Setting the timepoints from which the different models derive  
    timepoints_to_perform = rev(colnames(dataset_full)[as.numeric(colnames(dataset_full)) < end_timepoint])
    
    ###################### Leave one out Null #######################################################
    
    # Calculate mean accuracy of prediction on final timepoint, via Leave One Out method of splitting Train-Test sets
    LOO_multilogreg_Null <- function (dataset_full,timepoint,end_timepoint,meta_file,splitting_times){
      
      logreg <- as.data.frame(matrix(NA,nrow = nrow(dataset_full),ncol = 2))
      colnames(logreg)= c(paste('Cluster_at',end_timepoint,sep = '_'), timepoint)
      rownames(logreg)= rownames(dataset_full)
      
      for (modeling_sample in (rownames(dataset_full))){
        
        logreg[modeling_sample, 2] = ifelse(test = length(dataset_full[modeling_sample,timepoint]), yes = dataset_full[modeling_sample,timepoint], no = NA)
        logreg[modeling_sample, 1] = ifelse(test = length(dataset_full[modeling_sample,as.character(end_timepoint)]), yes = dataset_full[modeling_sample,as.character(end_timepoint)], no = NA)
        
      }
      
      logreg = logreg[complete.cases(logreg),]
      
      acc <- c()
      tra <- c()
      form = as.formula(paste(paste('Cluster_at',end_timepoint,sep = '_'),'~.',sep = ''))
      for (repeating_time in 1:splitting_times){
        for (j in 1:nrow(logreg)){
          trainset = logreg[1:nrow(logreg)!=j,]
          testset = logreg[j,]
          
          prediction_model <- multinom(formula = form , data = trainset ,censored = T, model = T)
          
          training_set_accuracies <- prediction_model %>% predict(trainset)
          test_set_predictions <- prediction_model %>% predict(testset)
          acc[j] = (sum(test_set_predictions == testset[,1])/ nrow(testset)) *100
          tra[j] = (sum(test_set_predictions == trainset[,1])/nrow(trainset))*100
          
        }
      }
      return (c(mean(tra),mean(acc)))
    }
    
    ###################### Stratified Train/Test splits Null ########################################
    
    # Calculate accuracy of prediction on final timepoint, as a mean of 100 accuracies with different
    # Train-Test stratified splits. Stratfication is performed on all metadata categories 
    multilogreg_stratified_Null <- function(dataset_full,optimal_splitting_percentage,timepoint,end_timepoint,meta_file, splitting_times){
      
      logreg <- as.data.frame(matrix(NA,nrow = nrow(dataset_full),ncol = 2))
      colnames(logreg)= c(paste('Cluster_at',end_timepoint,sep = '_'),timepoint)
      rownames(logreg)= rownames(dataset_full)
      
      for (modeling_sample in (rownames(dataset_full))){
        
        logreg[modeling_sample, 2] = ifelse(test = length(dataset_full[modeling_sample,as.character(timepoint)]), yes = dataset_full[modeling_sample,as.character(timepoint)], no = NA)
        logreg[modeling_sample, 1] = ifelse(test = length(dataset_full[modeling_sample,as.character(end_timepoint)]), yes = dataset_full[modeling_sample,as.character(end_timepoint)], no = NA)
        
      }
      logreg = logreg[complete.cases(logreg),]
      
      train_index <- createDataPartition(y = logreg[,paste('Cluster_at',end_timepoint,sep = '_')], p = optimal_splitting_percentage, list = F, times = splitting_times, groups = max(apply(X = logreg, MARGIN = 2, FUN = function(x){length(unique(x))})) )
      
      acc <- c()
      tra <- c()
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
    optimal_splitting_percentage <- function (dataset_full){
      all_calculated_accuracies <- c()
      best_splitting_percentage <- c()
      for (percentage_to_calculate in 65:90){
        percentage_to_calculate = percentage_to_calculate/100
        telika <- multilogreg_stratified_Null(dataset_full = dataset_full, optimal_splitting_percentage = percentage_to_calculate, timepoint = tail(timepoints_to_perform,1), end_timepoint = end_timepoint,meta_file = files,splitting_times = splitting_times)
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
    optimal_splitting_percentage <- optimal_splitting_percentage(dataset_full = dataset_full)
    
    ###################### Initializing variables Null ########################################
    Test_set_LOO_accuracy <- c()
    Train_set_LOO_accuracy <-c()
    Test_set_stratified_accuracy <- c()
    Train_set_stratified_accuracy <- c()
    trainsd <- c()
    testsd <- c()
    
    ###################### Calculate Accuracies on timepoints Null ############################
    
    
    
    # Loop to create a model from every timepoint to the final timepoint as selected from the first loop
    for (timepoint in timepoints_to_perform){
      
      LOOTeliko <- LOO_multilogreg_Null(dataset_full = dataset_full, timepoint = timepoint,end_timepoint = end_timepoint,meta_file = files, splitting_times = splitting_times)
      
      Train_set_LOO_accuracy = c(Train_set_LOO_accuracy,LOOTeliko[1])
      Test_set_LOO_accuracy = c(Test_set_LOO_accuracy,LOOTeliko[2])
      
      MultiLogReg_results <- multilogreg_stratified_Null(dataset_full = dataset_full, optimal_splitting_percentage = optimal_splitting_percentage,timepoint =  timepoint,meta_file =  files, end_timepoint = end_timepoint,splitting_times = splitting_times)
      
      
      Test_set_stratified_accuracy <- c(Test_set_stratified_accuracy,MultiLogReg_results[1])
      Train_set_stratified_accuracy <- c(Train_set_stratified_accuracy, MultiLogReg_results[2])
      trainsd <- c(trainsd, MultiLogReg_results[3])
      testsd <- c(testsd, MultiLogReg_results[4])
    }  
    
    
    ###################### Save accuracies on the lists ######################################
    
    # Adding to the lists in order to export the correct matrix
    Train_sets_LOO = rbind (Train_sets_LOO, c('Null',rep(0,(ncol(dataset_full)-length(Train_set_LOO_accuracy)-1)),Train_set_LOO_accuracy))
    Test_sets_LOO  = rbind (Test_sets_LOO,  c('Null',rep(0,(ncol(dataset_full)-length(Test_set_LOO_accuracy) -1)), Test_set_LOO_accuracy))
    Test_sets_stratified_split  = rbind (Test_sets_stratified_split,  c('Null',rep(0,(ncol(dataset_full)-length(Test_set_stratified_accuracy) -1)), Test_set_stratified_accuracy))
    Train_sets_stratified_split = rbind (Train_sets_stratified_split, c('Null',rep(0,(ncol(dataset_full)-length(Train_set_stratified_accuracy)-1)), Train_set_stratified_accuracy))
    
    
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
        LOO_multilogreg_One <- function (dataset_full,timepoint,end_timepoint,meta_file,Ncomb, splitting_times){
          
          logreg <- as.data.frame(matrix(NA,nrow = nrow(dataset_full),ncol = ncol(meta_file)))
          colnames(logreg)= c(paste('Cluster_at',end_timepoint,sep = '_'), timepoint , colnames(meta_file)[3:ncol(meta_file)])
          rownames(logreg)= rownames(dataset_full)
          
          for (modeling_sample in (rownames(dataset_full))){
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
            logreg[modeling_sample, 2] = ifelse(test = length(dataset_full[modeling_sample,timepoint]), yes = dataset_full[modeling_sample,timepoint], no = NA)
            logreg[modeling_sample, 1] = ifelse(test = length(dataset_full[modeling_sample,as.character(end_timepoint)]), yes = dataset_full[modeling_sample,as.character(end_timepoint)], no = NA)
          }
          
          logreg = logreg[complete.cases(logreg),]
          
          acc <- c()
          tra <- c()
          form = as.formula(paste(paste('Cluster_at',end_timepoint,sep = '_'),'~.',sep = ''))
          for (repeat_time in 1:splitting_times){
            for (j in 1:nrow(logreg)){
              trainset = logreg[1:nrow(logreg)!=j,]
              testset = logreg[j,]
              
              prediction_model <- multinom(formula = form , data = trainset ,censored = T, model = T)
              
              training_set_accuracies <- prediction_model %>% predict(trainset[,2:ncol(trainset)])
              test_set_predictions <- prediction_model %>% predict(testset[,2:ncol(testset)])
              acc[j] = (sum(test_set_predictions == testset[,1])/ nrow(testset)) *100
              tra[j] = (sum(test_set_predictions == trainset[,1])/nrow(trainset))*100
              
            }
          }        
          return (c(mean(tra),mean(acc)))
        }
        
        ###################### Stratified Train/Test splits Combinations ########################################
        
        # Calculate accuracy of prediction on final timepoint, as a mean of 100 accuracies with different
        # Train-Test stratified splits. Stratfication is performed on all metadata categories 
        multilogreg_stratified_One <- function(dataset_full,optimal_splitting_percentage,timepoint,end_timepoint,meta_file,Ncomb,splitting_times){
          
          logreg <- as.data.frame(matrix(NA,nrow = nrow(dataset_full),ncol = ncol(meta_file)))
          colnames(logreg)= c(paste('Cluster_at',end_timepoint,sep = '_'),timepoint,colnames(meta_file)[3:ncol(meta_file)])
          rownames(logreg)= rownames(dataset_full)
          
          for (modeling_sample in (rownames(dataset_full))){
            
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
            logreg[modeling_sample, 2] = ifelse(test = length(dataset_full[modeling_sample,timepoint]), yes = dataset_full[modeling_sample,timepoint], no = NA)
            logreg[modeling_sample, 1] = ifelse(test = length(dataset_full[modeling_sample,as.character(end_timepoint)]), yes = dataset_full[modeling_sample,as.character(end_timepoint)], no = NA)
            
            
          }
          logreg = logreg[complete.cases(logreg),]
          
          train_index <- createDataPartition(y = logreg[,paste('Cluster_at',end_timepoint,sep = '_')], p = optimal_splitting_percentage, list = F, times = splitting_times, groups = max(apply(X = logreg, MARGIN = 2, FUN = function(x){length(unique(x))})) )
          
          acc <- c()
          tra <- c()
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
        optimal_splitting_percentage <- function (dataset_full){
          all_calculated_accuracies <- c()
          best_splitting_percentage <- c()
          for (percentage_to_calculate in 65:90){
            percentage_to_calculate = percentage_to_calculate/100
            telika <- multilogreg_stratified_One(dataset_full = dataset_full, optimal_splitting_percentage = percentage_to_calculate, timepoint = tail(timepoints_to_perform,1), end_timepoint = end_timepoint,meta_file = files,Ncomb = Ncomb,splitting_times = splitting_times)
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
        optimal_splitting_percentage <- optimal_splitting_percentage(dataset_full = dataset_full)
        
        ###################### Initializing variables Combinations ########################################
        Test_set_LOO_accuracy <- c()
        Train_set_LOO_accuracy <-c()
        Test_set_stratified_accuracy <- c()
        Train_set_stratified_accuracy <- c()
        trainsd <- c()
        testsd <- c()
        
        ###################### Calculate Accuracies on timepoints Combinations ############################
        
        # Loop to create a model from every timepoint to the final timepoint as selected from the first loop
        for (timepoint in timepoints_to_perform){
          
          LOOTeliko <- LOO_multilogreg_One(dataset_full = dataset_full, timepoint = timepoint,end_timepoint = end_timepoint,meta_file = files,Ncomb = Ncomb, splitting_times = splitting_times)
          
          Train_set_LOO_accuracy = c(Train_set_LOO_accuracy,LOOTeliko[1])
          Test_set_LOO_accuracy = c(Test_set_LOO_accuracy,LOOTeliko[2])
          
          MultiLogReg_results <- multilogreg_stratified_One(dataset_full = dataset_full, optimal_splitting_percentage = optimal_splitting_percentage,timepoint =  timepoint,meta_file =  files, end_timepoint = end_timepoint, Ncomb = Ncomb,splitting_times = splitting_times)
          
          
          Test_set_stratified_accuracy <- c(Test_set_stratified_accuracy,MultiLogReg_results[1])
          Train_set_stratified_accuracy <- c(Train_set_stratified_accuracy, MultiLogReg_results[2])
          trainsd <- c(trainsd, MultiLogReg_results[3])
          testsd <- c(testsd, MultiLogReg_results[4])
        }  
        
        
        ###################### Save accuracies on the lists ######################################
        
        # Select the effects that created the model
        effect = paste(effector_combinations[combination,],collapse = ' & ')
        
        # Adding to the lists in order to export the correct matrix
        Train_sets_LOO = rbind (Train_sets_LOO, c(effect,rep(0,(ncol(dataset_full)-length(Train_set_LOO_accuracy)-1)),   Train_set_LOO_accuracy))
        Test_sets_LOO  = rbind (Test_sets_LOO,  c(effect,rep(0,(ncol(dataset_full)-length(Test_set_LOO_accuracy)-1)),    Test_set_LOO_accuracy))
        Test_sets_stratified_split  = rbind (Test_sets_stratified_split,  c(effect,rep(0,(ncol(dataset_full)-length(Test_set_stratified_accuracy)-1)), Test_set_stratified_accuracy))
        Train_sets_stratified_split = rbind (Train_sets_stratified_split, c(effect,rep(0,(ncol(dataset_full)-length(Train_set_stratified_accuracy)-1)),Train_set_stratified_accuracy))
        
      }
    }
    
    
    
    
    
    ###################### Leave one out All #######################################
    # Calculate mean accuracy of prediction on final timepoint, via Leave One Out method of splitting Train-Test sets
    LOO_multilogreg <- function (dataset_full,timepoint, end_timepoint, splitting_times){
      
      
      # Create the matrix with columns the state on the timepoint, the metadata and as a response the state on the last timepoint
      logreg <- as.data.frame(matrix(NA,nrow = nrow(dataset_full),ncol = ncol(meta_file)))
      colnames(logreg)= c(paste('Cluster_at',end_timepoint,sep = '_'), timepoint , colnames(meta_file)[3:ncol(meta_file)])
      rownames(logreg)= rownames(dataset_full)
      
      for (modeling_sample in (rownames(dataset_full))){
        subTab = meta_file[meta_file[,1] %in% modeling_sample & meta_file[,2] == as.numeric(end_timepoint), 3:ncol(meta_file)]
        if (nrow(subTab) == 1){
          logreg[modeling_sample, colnames(meta_file)[3:ncol(meta_file)]] = subTab[1,]
        }
        logreg[modeling_sample, 2] = ifelse(test = length(dataset_full[modeling_sample,timepoint]), yes = dataset_full[modeling_sample,timepoint], no = NA)
        logreg[modeling_sample, 1] = ifelse(test = length(dataset_full[modeling_sample,as.character(end_timepoint)]), yes = dataset_full[modeling_sample,as.character(end_timepoint)], no = NA)
      }
      
      # Remove NAs
      logreg = logreg[complete.cases(logreg),]
      
      # Set the response variable
      fo = as.formula(paste(paste('Cluster_at',end_timepoint,sep = '_'),'~.',sep = ''))
      
      acc <- c()
      tra <- c()
      for (repeat_time in 1:splitting_times){
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
      }
      
      return (c(mean(tra),mean(acc)))
    }
    
    ###################### Stratified Train/Test splits All ########################
    
    # Calculate accuracy of prediction on final timepoint, as a mean of 100 accuracies with different
    # Train-Test stratified splits. Stratfication is performed on all metadata categories 
    multilogreg_stratified <- function(dataset_full,optimal_splitting_percentage,timepoint,times, end_timepoint,splitting_times){
      
      # Create the matrix with columns the state on the timepoint, the metadata and as a response the state on the last timepoint
      logreg <- as.data.frame(matrix(NA,nrow = nrow(dataset_full),ncol = ncol(meta_file)))
      colnames(logreg)= c(paste('Cluster_at',end_timepoint,sep = '_'),timepoint,colnames(meta_file)[3:ncol(meta_file)])
      rownames(logreg)= rownames(dataset_full)
      
      for (modeling_sample in (rownames(dataset_full))){
        subTab = meta_file[meta_file[,1] %in% modeling_sample & meta_file[,2] == as.numeric(end_timepoint), 3:ncol(meta_file)]
        if (nrow(subTab) == 1){
          logreg[modeling_sample, colnames(meta_file)[3:ncol(meta_file)]] = subTab[1,]
        }
        logreg[modeling_sample, 2] = ifelse(test = length(dataset_full[modeling_sample,timepoint]), yes = dataset_full[modeling_sample,timepoint], no = NA)
        logreg[modeling_sample, 1] = ifelse(test = length(dataset_full[modeling_sample,as.character(end_timepoint)]), yes = dataset_full[modeling_sample,as.character(end_timepoint)], no = NA)
      }
      # Remove NAs
      logreg = logreg[complete.cases(logreg),]
      # Create different partitions to split train and test sets
      train_index <- createDataPartition(y = logreg[,paste('Cluster_at',end_timepoint,sep = '_')], p = optimal_splitting_percentage, list = F, times = splitting_times, groups = max(apply(X = logreg, MARGIN = 2, FUN = function(x){length(unique(x))})) )
      
      # Set the response variable
      fo = as.formula(paste(paste('Cluster_at',end_timepoint,sep = '_'),'~.',sep = ''))
      
      acc <- c()
      tra <- c()
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
    optimal_splitting_percentage <- function (dataset_full){
      all_calculated_accuracies <- c()
      best_splitting_percentage <- c()
      for (percentage_to_calculate in 65:90){
        percentage_to_calculate = percentage_to_calculate/100
        telika <- multilogreg_stratified(end_timepoint = end_timepoint,dataset_full = dataset_full, optimal_splitting_percentage = percentage_to_calculate, timepoint = rev(colnames(dataset_full))[2], splitting_times = splitting_times)
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
    optimal_splitting_percentage <- optimal_splitting_percentage(dataset_full = dataset_full)
    
    
    ###################### Initializing variables All #############################
    
    Test_set_LOO_accuracy <- c()
    Train_set_LOO_accuracy <-c()
    Test_set_stratified_accuracy <- c()
    Train_set_stratified_accuracy <- c()
    trainsd <- c()
    testsd <- c()
    
    
    ###################### Calculate Accuracies on timepoints All #################
    
    for (timepoint in timepoints_to_perform){
      
      
      LOOTeliko <- LOO_multilogreg(dataset_full = dataset_full, timepoint = timepoint,end_timepoint = end_timepoint, splitting_times = splitting_times)
      
      Train_set_LOO_accuracy = c(Train_set_LOO_accuracy,LOOTeliko[1])
      Test_set_LOO_accuracy = c(Test_set_LOO_accuracy,LOOTeliko[2])
      
      MultiLogReg_results <- multilogreg_stratified(dataset_full, optimal_splitting_percentage, timepoint, times = 100,end_timepoint = end_timepoint,splitting_times = splitting_times)
      
      
      Test_set_stratified_accuracy <- c(Test_set_stratified_accuracy,MultiLogReg_results[1])
      Train_set_stratified_accuracy <- c(Train_set_stratified_accuracy, MultiLogReg_results[2])
      trainsd <- c(trainsd, MultiLogReg_results[3])
      testsd <- c(testsd, MultiLogReg_results[4])
    }  
    
    ###################### Save accuracies on the lists ######################################
    Train_sets_LOO = rbind (Train_sets_LOO,  c('All',rep(0,(ncol(dataset_full)-length(Train_set_LOO_accuracy)-1)),Train_set_LOO_accuracy))
    Test_sets_LOO  = rbind (Test_sets_LOO,   c('All',rep(0,(ncol(dataset_full)-length(Test_set_LOO_accuracy)-1)),Test_set_LOO_accuracy))
    Test_sets_stratified_split  = rbind (Test_sets_stratified_split,   c('All',rep(0,(ncol(dataset_full)-length(Test_set_stratified_accuracy)-1)),Test_set_stratified_accuracy))
    Train_sets_stratified_split = rbind (Train_sets_stratified_split,  c('All',rep(0,(ncol(dataset_full)-length(Train_set_stratified_accuracy)-1)),Train_set_stratified_accuracy))
    Train_stratified_overtime_accuracy <- c(Train_set_stratified_accuracy,Train_stratified_overtime_accuracy)
    Test_stratified_overtime_accuracy  <- c(Test_set_stratified_accuracy ,Test_stratified_overtime_accuracy)
    Train_LOO_overtime_accuracy <- c(Train_set_LOO_accuracy,Train_LOO_overtime_accuracy)
    Test_LOO_overtime_accuracy  <- c(Test_set_LOO_accuracy, Test_LOO_overtime_accuracy)

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
  
  colnames(Train_sets_LOO) = c(paste('From timepoint ',rev(colnames(dataset_full)[1:ncol(dataset_full)-1]),sep = ''))
  colnames(Test_sets_LOO)  = c(paste('From timepoint ',rev(colnames(dataset_full)[1:ncol(dataset_full)-1]),sep = ''))
  colnames(Train_sets_stratified_split) = c(paste('From timepoint ',rev(colnames(dataset_full)[1:ncol(dataset_full)-1]),sep = ''))
  colnames(Test_sets_stratified_split)  = c(paste('From timepoint ',rev(colnames(dataset_full)[1:ncol(dataset_full)-1]),sep = ''))
  write.table(x = Train_sets_LOO, file = paste(output_dir,'All Accuracies of Training Sets LOO.tab', sep = '/') , sep = "\t",col.names =NA, row.names = TRUE,quote = FALSE)
  write.table(x = Test_sets_LOO,  file = paste(output_dir,'All Accuracies of Test Sets LOO.tab',     sep = '/') , sep = "\t",col.names =NA, row.names = TRUE,quote = FALSE)
  write.table(x = Train_sets_stratified_split, file = paste(output_dir,'All Accuracies of Training Sets Splitted.tab', sep = '/') , sep = "\t",col.names =NA, row.names = TRUE,quote = FALSE)
  write.table(x = Test_sets_stratified_split,  file = paste(output_dir,'All Accuracies of Test Sets Splitted.tab', sep = '/')  , sep = "\t",col.names =NA, row.names = TRUE,quote = FALSE)
  
  
  ###################### Barplots of Accuracies ############################################
  
  # Create a directory wherre the acuuracy related files will be placed
  dir.create(paste(output_dir,"Accuracies",sep = '/'),showWarnings = F)
  
  
  # Counting indices
  row <- 1 ; coloumn <- 0
  
  
  # Calculate the maximum accuracy for every timepoint
  for (x in 1:ncol(Train_sets_LOO)){
    
    # Vector with the maximum accuracies and the metadata
    max_Train_sets_LOO <- c()
    meta_Train_sets_LOO <- c()
    max_Test_sets_LOO <- c()
    meta_Test_sets_LOO <- c()
    max_Train_sets_stratified_split <- c()
    meta_Train_sets_stratified_split <- c()
    max_Test_sets_stratified_split <- c()
    meta_Test_sets_stratified_split <- c()
    
    for (y in 1:ncol(Train_sets_LOO)){
      if ((coloumn+y) <= ncol(Train_sets_LOO)){
        max_Train_sets_LOO <- c( max_Train_sets_LOO, max(as.numeric(Train_sets_LOO[row:Npredictions,coloumn+y])))
        meta_Train_sets_LOO <- c(meta_Train_sets_LOO, names(which(Train_sets_LOO[row:Npredictions,coloumn+y]== as.character(max(as.numeric(Train_sets_LOO[row:Npredictions,coloumn+y]))))[1]))
        max_Test_sets_LOO <- c( max_Test_sets_LOO,max(as.numeric(Test_sets_LOO[row:Npredictions,coloumn+y])))
        meta_Test_sets_LOO <- c(meta_Test_sets_LOO, names(which(Test_sets_LOO[row:Npredictions,coloumn+y]== as.character(max(as.numeric(Test_sets_LOO[row:Npredictions,coloumn+y]))))[1]))
        max_Train_sets_stratified_split <- c( max_Train_sets_stratified_split, max(as.numeric(Train_sets_stratified_split[row:Npredictions,coloumn+y])))
        meta_Train_sets_stratified_split<- c(meta_Train_sets_stratified_split, names(which(Train_sets_stratified_split[row:Npredictions,coloumn+y]== as.character(max(as.numeric(Train_sets_stratified_split[row:Npredictions,coloumn+y]))))[1]))
        max_Test_sets_stratified_split <- c( max_Test_sets_stratified_split, max(as.numeric(Test_sets_stratified_split[row:Npredictions,coloumn+y])))
        meta_Test_sets_stratified_split <- c(meta_Test_sets_stratified_split, names(which(Test_sets_stratified_split[row:Npredictions,coloumn+y]== as.character( max(as.numeric(Test_sets_stratified_split[row:Npredictions,coloumn+y]))))[1]))
      }
    }
    
    # Save the accuracies in the heatmap matrix
    heatmap_matrix_statified[ncol(heatmap_matrix_LOO) - x +1,c(1:length(max_Test_sets_stratified_split))] <- max_Test_sets_stratified_split
    
    heatmap_matrix_LOO[ncol(heatmap_matrix_LOO) - x +1,c(1:length(max_Test_sets_LOO))] <- max_Test_sets_LOO
    
    # Formation of the metadata matrix 
    temp_df_meta <- data.frame(rbind(meta_Test_sets_LOO,meta_Test_sets_stratified_split, meta_Train_sets_LOO,meta_Train_sets_stratified_split))
    temp_df_meta <- data.frame(c("Test_sets_LOO","Test_sets_stratified_split","Train_sets_LOO","Train_sets_stratified_split"),temp_df_meta)
    colnames(temp_df_meta) <- c("rows",colnames(dataset_full_on_clusters)[1:(ncol(Train_sets_LOO)-coloumn)])
    
    # Unpivot temp_df_meta from wide to long format using the row names as the identifier
    temp_df_meta <-  data.frame(melt(temp_df_meta, cols = 2:ncol(temp_df_meta),id.vars ="rows"))
    colnames(temp_df_meta) <- c("Clusters","Timepoints","Percentages")
    
    # Convert Timepoints and Cluster columns into factors
    temp_df_meta[,"Clusters"] <- factor(temp_df_meta[,"Clusters"],levels = unique(temp_df_meta[,"Clusters"]))
    temp_df_meta[,"Timepoints"] <- factor( temp_df_meta[,"Timepoints"], levels=colnames(dataset_full_on_clusters)[1:(ncol(Train_sets_LOO)-coloumn)])
    
    
    # Formation of the accuracies matrix 
    temp_df <- data.frame(rbind(max_Test_sets_LOO,max_Test_sets_stratified_split, max_Train_sets_LOO,max_Train_sets_stratified_split))
    temp_df <- data.frame(c("Test_sets_LOO","Test_sets_stratified_split","Train_sets_LOO","Train_sets_stratified_split"),temp_df)
    colnames(temp_df) <- c("rows",colnames(dataset_full_on_clusters)[1:(ncol(Train_sets_LOO)-coloumn)])
    
    # Unpivot temp_df from wide to long format using the row names as the identifier
    temp_df <-  data.frame(melt(temp_df, cols = 2:ncol(temp_df),id.vars ="rows"))
    colnames(temp_df) <- c("Clusters","Timepoints","Percentages")
    
    # Convert Timepoints and Cluster columns into factors
    temp_df[,"Clusters"] <- factor(temp_df[,"Clusters"],levels = unique(temp_df[,"Clusters"]))
    temp_df[,"Timepoints"] <- factor( temp_df[,"Timepoints"], levels=colnames(dataset_full_on_clusters)[1:(ncol(Train_sets_LOO)-coloumn)])
    temp_df <- data.frame(temp_df, meta=as.factor(temp_df_meta[,3]))
    
    # Create a vector with the colours for the barplot
    color <-colours_ploting[1:nlevels(as.factor(temp_df[,"Clusters"]))]
    
    # Barplot of the Accuracies
    barplot <-  ggplot(data=temp_df, aes(x=Timepoints, y=Percentages)) +
      scale_y_continuous(limits=c(0,100),breaks = c(0,20,40,60,80,100))+
      geom_bar(aes(fill= Clusters),stat="identity", position="dodge",alpha=0.7)+
      geom_text(aes(x=Timepoints, y=Percentages,label=meta, group = Clusters),position = position_dodge(width = .9),angle=90,hjust=2.5)+
      ylab("Acurracy")+
      ggtitle(paste("Accuracy achieved on Timepoint:",rev(colnames(dataset_full_on_clusters))[x]))+ 
      scale_fill_manual(breaks=c(levels(factor(temp_df[,"Clusters"],levels = unique(temp_df[,"Clusters"])))), values=color)+
      guides(fill=guide_legend(title=""))+ 
      geom_hline(yintercept = (1/nlevels(as.factor(dataset_full_on_clusters[,rev(colnames(dataset_full_on_clusters))[x]])))*100, linetype = "dashed")+
      theme_classic()+
      theme(legend.title.align = 0.5,plot.title = element_text(hjust = 0.5))+
      scale_fill_manual(values = c('paleturquoise4','paleturquoise3', 'darkolivegreen2', 'darkolivegreen3'),guide="none")
    
    # Calculate the initial line for the next iteration
    row <- row+Npredictions
    # Calculate the initial column for the next iteration
    coloumn=coloumn+1
    
    # Print the barplots
    ggsave(paste(paste(output_dir,"Accuracies",sep = '/'),paste(paste('Accuracies of Models on Timepoint', rev(colnames(dataset_full_on_clusters))[coloumn],sep = ' '),'pdf', sep='.'),sep = '/'),barplot)
    jpeg(filename = paste(paste(output_dir,"Accuracies",sep = '/'),paste(paste('Accuracies of Models on Timepoint', rev(colnames(dataset_full_on_clusters))[coloumn],sep = ' '),'jpeg', sep='.'),sep = '/'))
    print(barplot)
    dev.off()
    
  }
  
  
  # Function for the Heatmaps
  heatmap <- function (heatmap_matrix,category) {
    # Name the rows and columns of the heatmap matrix
    rownames(heatmap_matrix) <- colnames(dataset_full_on_clusters)
    colnames(heatmap_matrix) <- colnames(dataset_full_on_clusters)
    
    # Unpivot heatmap_matrix from wide to long format using the rownames as the identifier
    heatmap_matrix <- t(heatmap_matrix)
    heatmap_matrix <- data.frame(rownames(heatmap_matrix),heatmap_matrix)
    colnames(heatmap_matrix)[2:ncol(heatmap_matrix)] <- colnames(dataset_full_on_clusters)
    melted_heatmap <- melt(heatmap_matrix, na.rm = TRUE)
    colnames(melted_heatmap) <- c("X","Y","Value")
    
    # Convert X and Y columns into factors
    melted_heatmap[,"X"] <- factor( melted_heatmap[,"X"] , levels=colnames(dataset_full_on_clusters))
    melted_heatmap[,"Y"] <- factor( melted_heatmap[,"Y"] , levels=colnames(dataset_full_on_clusters))
    
    # Heatmap
    heatmap_plot <-  ggplot(data = melted_heatmap, aes(X, Y, fill = Value))+
      geom_tile(color = "white")+
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 50, limit = c(0,100), space = "Lab", name="Accuracy") +
      theme_minimal()+
      ylab("Timepoints")+
      xlab("Timepoints")+
      coord_fixed()+
      geom_text(aes(X, Y, label = round(Value,2)), color = "black", size = 4) +
      theme(plot.title = element_text(hjust = 0.5,size=20, face = "bold"),
            legend.box.just = "left",
            axis.text.x = element_text(size = 12, face = "bold"),
            axis.text.y = element_text(size = 12, face = "bold"),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            legend.justification = c(0, 1),
            legend.position = c(0.7, 0.3),
            legend.direction = "horizontal")+
      guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))+
      ggtitle(paste0("Heatmap of Accuracies (",category,")"))
    
    # Print the heatmap
    ggsave(paste(paste(output_dir,"Accuracies",sep = '/'),paste0('Heatmap(',category,').pdf'),sep = '/'),heatmap_plot)
    jpeg(filename = paste(paste(output_dir,"Accuracies",sep = '/'),paste0('Heatmap(',category,').jpeg'),sep = '/'))
    print(heatmap_plot)
    dev.off()
  }
  
  # Print the LOO heatmap 
  heatmap(heatmap_matrix_LOO,"LOO")
  
  # Print the stratified Split heatmap 
  heatmap(heatmap_matrix_statified,"Stratified Split")
  
  ###################### Calculate random estimators performance ##########################################
  
  random_estimator = 100/rev(apply (X = dataset_full,MARGIN = 2,FUN = function(x){max(x,na.rm = T)}))
  
  ###################### Write files of Best Accuracies Test_sets_LOO ###########################################
  
  TotimepointIndeces = c(0:ncol(Test_sets_LOO))*Npredictions
  maxaccuracies <- c()
  metadata <- c()
  for (i in 1:(length(TotimepointIndeces)-1)){
  
    maxaccuracies[i] <-  max(as.numeric(Test_sets_LOO[TotimepointIndeces[i]:TotimepointIndeces[i+1],i]))
    metadata[i] <-  names(which(Test_sets_LOO[TotimepointIndeces[i]:TotimepointIndeces[i+1],i]== as.character(max(as.numeric(Test_sets_LOO[TotimepointIndeces[i]:TotimepointIndeces[i+1],i]))))[1])
    
     
    #maxaccuracies[i]   = max(Test_sets_LOO[TotimepointIndeces[i]:TotimepointIndeces[i+1],i])
    #metadata[i] = rownames(Test_sets_LOO)[which.max(Test_sets_LOO[TotimepointIndeces[i]:TotimepointIndeces[i+1],i])]
    
  }
  write.table(x = rbind(paste(paste('Maximum Accuracy for Timepoint', rev(colnames(dataset_full)[2:ncol(dataset_full)])), c(rev(colnames(dataset_full))[2:ncol(dataset_full)]), sep = ' from Timepoint '),maxaccuracies,metadata,random_estimator[1:length(random_estimator)-1]),
              file = paste(output_dir,'Maximum_Accuracies_of_Test_sets_LOO.tab',sep = '/'), sep = "\t",col.names =NA, row.names = TRUE,quote = FALSE)
  

  ###################### Write files of Best Accuracies TestSplits ########################################
  
  maxaccuracies <- c()
  metadata <- c()
  for (i in 1:(length(TotimepointIndeces)-1)){
    
    maxaccuracies[i] <-  max(as.numeric(Test_sets_stratified_split[TotimepointIndeces[i]:TotimepointIndeces[i+1],i]))
    metadata[i] <-  names(which(Test_sets_stratified_split[TotimepointIndeces[i]:TotimepointIndeces[i+1],i]== as.character(max(as.numeric(Test_sets_stratified_split[TotimepointIndeces[i]:TotimepointIndeces[i+1],i]))))[1])
    
    # maxaccuracies[i]   = max(Test_sets_stratified_split[TotimepointIndeces[i]:TotimepointIndeces[i+1],i])
    # metadata[i] = rownames(Test_sets_stratified_split)[which.max(Test_sets_stratified_split[TotimepointIndeces[i]:TotimepointIndeces[i+1],i])]
    
  }
  write.table(x = rbind(paste(paste('Maximum Accuracy for Timepoint', rev(colnames(dataset_full)[2:ncol(dataset_full)])), c(rev(colnames(dataset_full))[2:ncol(dataset_full)]), sep = ' from Timepoint '),maxaccuracies,metadata,random_estimator[1:length(random_estimator)-1])
              ,file = paste(output_dir,'Maximum_Accuracies_of_Stratified_Tests.tab',sep = '/'), sep = "\t",col.names =NA, row.names = TRUE,quote = FALSE)
  
  
  ###################### Write files of Best Accuracies Train_sets_LOO ########################################
  
  maxaccuracies <- c()
  metadata <- c()
  for (i in 1:(length(TotimepointIndeces)-1)){
    
    maxaccuracies[i] <-  max(as.numeric(Train_sets_LOO[TotimepointIndeces[i]:TotimepointIndeces[i+1],i]))
    metadata[i] <-  names(which(Train_sets_LOO[TotimepointIndeces[i]:TotimepointIndeces[i+1],i]== as.character(max(as.numeric(Train_sets_LOO[TotimepointIndeces[i]:TotimepointIndeces[i+1],i]))))[1])
    
    # maxaccuracies[i]   = max(Train_sets_LOO[TotimepointIndeces[i]:TotimepointIndeces[i+1],i])
    # metadata[i] = rownames(Train_sets_LOO)[which.max(Train_sets_LOO[TotimepointIndeces[i]:TotimepointIndeces[i+1],i])]
    
  }
  write.table(x = rbind(paste(paste('Maximum Accuracy for Timepoint', rev(colnames(dataset_full)[2:ncol(dataset_full)])), c(rev(colnames(dataset_full))[2:ncol(dataset_full)]), sep = ' from Timepoint '),maxaccuracies,metadata,random_estimator[1:length(random_estimator)-1])
              ,file = paste(output_dir,'Maximum_Accuracies_of_Train_sets_LOO.tab',sep = '/'), sep = "\t",col.names =NA, row.names = TRUE,quote = FALSE)
  

  ###################### Write files of Best Accuracies TrainSplits ########################################
  
  maxaccuracies <- c()
  metadata <- c()
  
  for (i in 1:(length(TotimepointIndeces)-1)){
    
    maxaccuracies[i] <-  max(as.numeric(Train_sets_stratified_split[TotimepointIndeces[i]:TotimepointIndeces[i+1],i]))
    metadata[i] <-  names(which(Train_sets_stratified_split[TotimepointIndeces[i]:TotimepointIndeces[i+1],i]== as.character(max(as.numeric(Train_sets_stratified_split[TotimepointIndeces[i]:TotimepointIndeces[i+1],i]))))[1])
    
    # maxaccuracies[i]   = max(Train_sets_stratified_split[TotimepointIndeces[i]:TotimepointIndeces[i+1],i])
    # metadata[i] = rownames(Train_sets_stratified_split)[which.max(Train_sets_stratified_split[TotimepointIndeces[i]:TotimepointIndeces[i+1],i])]
    
  }
  
  write.table(x = rbind(paste(paste('Maximum Accuracy for Timepoint', rev(colnames(dataset_full)[2:ncol(dataset_full)])), c(rev(colnames(dataset_full))[2:ncol(dataset_full)]), sep = ' from Timepoint '),maxaccuracies,metadata,random_estimator[1:length(random_estimator)-1]),
              file = paste(output_dir,'Maximum_Accuracies_of_TrainSplits.tab',sep = '/'), sep = "\t",col.names =NA, row.names = TRUE,quote = FALSE)
  
  ###################### Perform Chi-square analysis to check metadata influence on T0 #####################
  
  # Create a matrix to check for the effect of the metadata on the first timepoint (T0)
  chimatrix = matrix (NA, ncol = length(effectors)+1, nrow = nrow(dataset_full))
  
  # Assign the rownames to the matrix
  rownames(chimatrix) = rownames(dataset_full)
  
  # Assign the values on the matrix
  chimatrix[,1] = dataset_full[,1]
  
  for (chirow in rownames(dataset_full)){
    
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
  
  
  
  print (paste('The transitions among timepoints ',ifelse(test = Markovian_property, yes = 'are Markovian', no = 'are NOT Markovian'),sep = ' '))
  print (paste('For the clustering fate on the first timepoint', paste(ifelse(test = Significant_metadata, yes = Significant_metadata, no = 'No'), 'metadata are significant.', sep = ' '), sep = ' '))
  print (' Analysis Completed ')
  print ('_____________________________________________')
  print ('Results are saved in the preselected folders ')
  print ('_____________________________________________')
  
  
  
} else {
  print (paste('Cronos has already run with the exact same parameters. The output files are stored in',sub(pattern = './',replacement = '',x = directory),sep = ' '))
}
