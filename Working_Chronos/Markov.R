#################### NO CHANGE IS NEEDED HERE ###########################################

#### You can select everything (Ctrl+A) and press run (Ctrl+Enter)  

source('Clustering.R')

##########################################################################################
####################### MARKOVIAN CHAIN CHECK ############################################
##########################################################################################


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


infants<- samples_on_clusters[is.na(samples_on_clusters[,adult_timepoint_name]),1:ncol(samples_on_clusters)-1]
kati_allo <- matrix (NA, nrow = nrow(infants), ncol = 5)
colnames(kati_allo) = sort(as.numeric(colnames(infants)[colnames(infants)!='ot']))
rownames(kati_allo) = rownames(infants)
kati_allo
colnames(infants)
colnames(kati_allo)
for (row in rownames(infants)){
  kati_allo[row,1:2] = infants[row,1:2]
  kati_allo[row,3] = infants[row,6]
  kati_allo[row,4] = infants[row,3]
  kati_allo[row,5]= infants[row,4]
}
kati_allo

infants_on_clusters <- samples_on_clusters[is.na(samples_on_clusters[,adult_timepoint_name]),1:ncol(samples_on_clusters)-1]
independed_clusters <- cumsum(as.vector(apply(X = infants,MARGIN = 2,FUN = function(x){max(x,na.rm = T)})))
for (i in 2:ncol(kati_allo)){
  for (j in 1:max(kati_allo[,i],na.rm = T)){
    kati_allo[kati_allo[,i]==j,i] <- independed_clusters[i-1]+j
  }
}

kati_allo
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
  girna <- matrix(NA, nrow = 3,ncol = length(dof))
  rownames(girna) = c('Markov_Property_per_transition','X2', 'Degrees_of_freedom')
  girna [1,] <- Markovian_Property_Per_Timepoints
  girna [2,] <- Ss
  girna [3,] <- dof
  return (girna)
}

Markovian_Property <- Markov_Test_X2(infants = infants)

Markovian_Property_renamed <-Markov_Test_X2(infants = infants_on_clusters)


######################### END OF SECTION #################################################

##########################################################################################
####################### CLAIM TRANSITION MATRIX ##########################################
##########################################################################################

transition_matrix <- matrix(0, nrow = max(kati_allo, na.rm = T), ncol = max(kati_allo,na.rm = T))

for (i in 1:(ncol(kati_allo)-1)){
  for (j in min(kati_allo[,i],na.rm = T):max(kati_allo[,i],na.rm = T)){
    for (k in min(kati_allo[,i+1],na.rm = T):max(kati_allo[,i+1],na.rm = T)){
  
        transition_matrix[j,k] = sum(kati_allo[kati_allo[,i]==j,(i+1)]==k,na.rm = T)/length(kati_allo[kati_allo[,i]==j,(i+1)])
      }
    }
  }

#transition_matrix


##########################################################################################

##########################################################################################
######################## SPECIFY TAXA ON CLUSTERS   ######################################
##########################################################################################


taxa_per_cluster <- function(taxa_matrix,samples_on_clusters,timepoint_list, representation_method){
  taxa_clusters <- list()
  for (i in colnames(samples_on_clusters)){
    for (j in min(samples_on_clusters[,i], na.rm = T):max(samples_on_clusters[,i],na.rm = T)){
      timiclj = samples_on_clusters[!is.na(samples_on_clusters[,i]),i]
      metatimepoint= meta_file[meta_file[,"Timepoint"]==i,]
      
      if (representation_method == 'median'){
        taxa_clusters[[paste(as.character(i),as.character(j),sep = ' cluster ')]]= as.matrix(apply(X = (taxa_matrix[rownames(metatimepoint[match(metatimepoint[,'Sample'],(as.vector(names(timiclj[timiclj==j]))),nomatch = F),]),]),MARGIN = 2,FUN = median))
      }
      
      else if (representation_method == 'mean'){
        taxa_clusters[[paste(as.character(i),as.character(j),sep = ' cluster ')]]= as.matrix(apply(X = (taxa_matrix[rownames(metatimepoint[match(metatimepoint[,'Sample'],(as.vector(names(timiclj[timiclj==j]))),nomatch = F),]),]),MARGIN = 2,FUN = mean))
      }
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

taxa_clusters <- taxa_per_cluster(taxa_matrix = taxa_matrix, samples_on_clusters = samples_on_clusters, timepoint_list = timepoint_list, representation_method = representation_method)

samples_on_clusters = samples_on_clusters[,names(timepoint_list)[names(timepoint_list)!='ot']]
samples_on_clusters = samples_on_clusters[,c(as.character(sort(as.numeric(colnames(samples_on_clusters)[colnames(samples_on_clusters)!=adult_timepoint_name]),decreasing = F,na.last = T)),adult_timepoint_name)]


infants<- samples_on_clusters[is.na(samples_on_clusters[,adult_timepoint_name]),1:ncol(samples_on_clusters)-1]

######################### END OF SECTION #################################################

################ WRITE TAB DELIMITED FILES WITH THE OUTPUTS  #############################
##########################################################################################

dir.create(dir_with_files, showWarnings = F)


colnames(taxa_clusters)<-lapply(X = colnames(taxa_clusters), FUN = function(x){paste('Timepoint',x,sep = ' ')})

write.csv(x = samples_on_clusters, file = paste(dir_with_files, "Samples_in_Timepoint-specific_Clusters.csv",sep = '/'), row.names = T)
write.csv(x = taxa_clusters,       file = paste(dir_with_files, "Taxonomic_profile_of_clusters.csv", sep = '/'), row.names = T)

######################### END OF SECTION #################################################


############## PRACTICE ############################

setwd('/home/arislitos/Working_Chronos/Chronos_almost/')
arnitiko <- read.csv(file = 'arnitiko_test.csv', sep = ',', header = T, check.names = F)
thetiko2 <- read.csv(file = 'thetiko2.csv', sep = ',', header = T, check.names = F)
plin <- Markov_Test_X2(arnitiko)[1,]
diko_mou <- Markov_Test_X2(infants = infants)[1,]
sin <- Markov_Test_X2(thetiko2)[1,]
l1= c('Positive Control')
l2 = c('Negative Control')

############## COMMENTS ############################


