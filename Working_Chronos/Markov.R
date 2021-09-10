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

infants_on_clusters <- samples_on_clusters[is.na(samples_on_clusters[,adult_timepoint_name]),1:ncol(samples_on_clusters)-1]
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
    Markovian_Property_Per_Timepoints[i] = Ss[i] < qchisq(.95,df = dof[i])
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

all(Markovian_Property['X2',]==Markovian_Property_renamed['X2',])

######################### END OF SECTION #################################################

##########################################################################################
####################### CLAIM TRANSITION MATRIX ##########################################
##########################################################################################

transition_matrix <- matrix(0, nrow = max(infants_on_clusters, na.rm = T), ncol = max(infants_on_clusters,na.rm = T))

for (i in 1:(ncol(infants_on_clusters)-1)){
  for (j in min(infants_on_clusters[,i],na.rm = T):max(infants_on_clusters[,i],na.rm = T)){
    for (k in min(infants_on_clusters[,i+1],na.rm = T):max(infants_on_clusters[,i+1],na.rm = T)){
  
        transition_matrix[j,k] = sum(infants_on_clusters[infants_on_clusters[,i]==j,(i+1)]==k,na.rm = T)/length(infants_on_clusters[infants_on_clusters[,i]==j,(i+1)])
      }
    }
  }

transition_matrix



############## PRACTICE ############################



############## COMMENTS ############################