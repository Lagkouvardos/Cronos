setwd('Working_Chronos/Chronos_almost/')

ptm <- proc.time()
source('Markov.R')

###################### Initializing Lists to save metrics of modeling ############################
# Initializing lists to write accuracies of models
trainLOO <- c()
testLOO  <- c()
trainStr <- c()
testStr  <- c()

# Naming the effectors to be imputed on the model
effectors = colnames(meta_file)[3:ncol(meta_file)]

###################### Perform Tasks on all Timepoints ###############################

# Loop to create a model for each timepoint on the dataset
for (end_timepoint in head(as.numeric(rev(colnames(infants))),n = (ncol(infants))-1)){
  
  
# Setting the timepoints from which the different models derive  
timepoints_to_perform = rev(colnames(infants)[as.numeric(colnames(infants)) < end_timepoint])

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
        
        for (deigma in (rownames(infants))){
          subTab = meta_file[meta_file[,1] %in% deigma & meta_file[,2] == as.numeric(end_timepoint), 3:ncol(meta_file)]
          if (Ncomb > 1 ){
            if(nrow(subTab) == 1){
              logreg[deigma, colnames(meta_file)[3:ncol(meta_file)]] = subTab
            }
          }
          else if (Ncomb == 1){
            if (length(subTab) == 1){
              logreg[deigma, colnames(meta_file)[3:ncol(meta_file)]] = subTab
            }
          }
          logreg[deigma, 2] = ifelse(test = length(infants[deigma,timepoint]), yes = infants[deigma,timepoint], no = NA)
          logreg[deigma, 1] = ifelse(test = length(infants[deigma,as.character(end_timepoint)]), yes = infants[deigma,as.character(end_timepoint)], no = NA)
        }
        
        logreg = logreg[complete.cases(logreg),]
        
        acc <- c()
        tra <- c()
        form = as.formula(paste(paste('Cluster_at',end_timepoint,sep = '_'),'~.',sep = ''))
        for (j in 1:nrow(logreg)){
          trainset = logreg[1:nrow(logreg)!=j,]
          testset = logreg[j,]
          
          apotelesmata <- multinom(formula = form , data = trainset ,censored = T, model = T)
          
          giatotrain <- apotelesmata %>% predict(trainset[,2:ncol(trainset)])
          provlepseis <- apotelesmata %>% predict(testset[,2:ncol(testset)])
          acc[j] = (sum(provlepseis == testset[,1])/ nrow(testset)) *100
          tra[j] = (sum(provlepseis == trainset[,1])/nrow(trainset))*100
          
        }
        
        return (c(mean(tra),mean(acc)))
      }
      
      ###################### Stratified Train/Test splits Combinations ########################################
      
      # Calculate accuracy of prediction on final timepoint, as a mean of 100 accuracies with different
      # Train-Test stratified splits. Stratfication is performed on all metadata categories 
      multilogreg_stratified_One <- function(infants,pososto,timepoint,end_timepoint,meta_file,Ncomb){
        
        logreg <- as.data.frame(matrix(NA,nrow = nrow(infants),ncol = ncol(meta_file)))
        colnames(logreg)= c(paste('Cluster_at',end_timepoint,sep = '_'),timepoint,colnames(meta_file)[3:ncol(meta_file)])
        rownames(logreg)= rownames(infants)
        
        for (deigma in (rownames(infants))){
          
          subTab = meta_file[meta_file[,1] %in% deigma & meta_file[,2] == as.numeric(end_timepoint), 3:ncol(meta_file)]
          
          if (Ncomb > 1){
            if (nrow(subTab) == 1){  
              logreg[deigma, colnames(meta_file)[3:ncol(meta_file)]] = subTab
            }
          }
          else if (Ncomb == 1){
            if (length(subTab) == 1){
              logreg[deigma, colnames(meta_file)[3:ncol(meta_file)]] = subTab
            }
          }
          logreg[deigma, 2] = ifelse(test = length(infants[deigma,timepoint]), yes = infants[deigma,timepoint], no = NA)
          logreg[deigma, 1] = ifelse(test = length(infants[deigma,as.character(end_timepoint)]), yes = infants[deigma,as.character(end_timepoint)], no = NA)
          
          
        }
        logreg = logreg[complete.cases(logreg),]
        
        train_index <- createDataPartition(y = logreg[,paste('Cluster_at',end_timepoint,sep = '_')], p = pososto, list = F, times = 100, groups = max(apply(X = logreg, MARGIN = 2, FUN = function(x){length(unique(x))})) )
        
        acc <- c()
        tra <- c()
        paragontes <- list()
        fo = as.formula(paste(paste('Cluster_at',end_timepoint,sep = '_'),'~.',sep = ''))
        for (j in 1:ncol(train_index)){
          trainset = logreg[train_index[,j],]
          testset = logreg[-train_index[,j],]
          if (all(apply(X = logreg, MARGIN = 2, FUN = function(x){length(unique(x))})== apply(X = trainset, MARGIN = 2, FUN = function(x){length(unique(x))}))){
            
            apotelesmata <- multinom(formula = fo, data = trainset ,censored = T, model = T)
            
            giatotrain <- apotelesmata %>% predict(trainset[,2:ncol(trainset)])
            provlepseis <- apotelesmata %>% predict(testset[,2:ncol(testset)])
            acc[j] = (sum(provlepseis == testset[,1])/ nrow(testset)) *100
            tra[j] = (sum(giatotrain == trainset[,1])/nrow(trainset))*100
          }
        }
        
        acc = na.omit(acc)
        tra = na.omit(tra)
        
        accteststr <- mean(acc)
        acctrainstr = mean(tra)
        trainsd <- sd(tra)
        testsd <- sd(acc)
        
        return (c(accteststr,acctrainstr, trainsd,testsd))
      }
      
      ###################### Calculate optimal split percentage ############################################
      # Calculate the optimal split percentage between .65 and .90 for predicting the final timepoint on the
      # exact previous one to be used in all splits for all timepoints. The criterion is the best accuracy  
      # as derived by the previous function with argument the different percentages.
      bestpososto <- function (infants){
        olataaccuracies <- c()
        sepoiopososto <- c()
        for (poso in 65:90){
          poso = poso/100
          telika <- multilogreg_stratified_One(infants = infants, pososto = poso, timepoint = tail(timepoints_to_perform,1), end_timepoint = end_timepoint,meta_file = files,Ncomb = Ncomb)
          accteststr = telika[1]
          acctrainstr = telika[2]
          trainsd = telika[3]
          testsd = telika[4]
          olataaccuracies= c(olataaccuracies,accteststr)
          sepoiopososto = c(sepoiopososto,poso)
        }
        
        return (sepoiopososto[which.max(olataaccuracies)])
      }
      
      # Assign the best percentage on a variable to use on the final calculations
      to_kalitero_pososto <- bestpososto(infants = infants)
      
      ###################### Initializing variables Combinations ########################################
      acctest <- c()
      acctrain <-c()
      accteststr <- c()
      acctrainstr <- c()
      trainsd <- c()
      testsd <- c()
      paragontes <- list()
      
      ###################### Calculate Accuracies on timepoints Combinations ############################
      
      # Loop to create a model from every timepoint to the final timepoint as selected from the first loop
      for (timepoint in timepoints_to_perform){
        
        LOOTeliko <- LOO_multilogreg_One(infants = infants, timepoint = timepoint,end_timepoint = end_timepoint,meta_file = files,Ncomb = Ncomb)
        
        acctrain = c(acctrain,LOOTeliko[1])
        acctest = c(acctest,LOOTeliko[2])
        
        Splitteliko <- multilogreg_stratified_One(infants = infants, pososto = to_kalitero_pososto,timepoint =  timepoint,meta_file =  files, end_timepoint = end_timepoint, Ncomb = Ncomb)
        
        
        accteststr <- c(accteststr,Splitteliko[1])
        acctrainstr <- c(acctrainstr, Splitteliko[2])
        trainsd <- c(trainsd, Splitteliko[3])
        testsd <- c(testsd, Splitteliko[4])
      }  
      
      
      ###################### Save accuracies on the lists ######################################
      
      # Select the effects that created the model
      effect = paste(effector_combinations[combination,],collapse = ' & ')
      
      # Adding to the lists in order to export the correct matrix
      trainLOO = rbind (trainLOO, c(effect,rep(0,(ncol(infants)-length(acctrain)-1)),   acctrain))
      testLOO  = rbind (testLOO,  c(effect,rep(0,(ncol(infants)-length(acctest)-1)),    acctest))
      testStr  = rbind (testStr,  c(effect,rep(0,(ncol(infants)-length(accteststr)-1)), accteststr))
      trainStr = rbind (trainStr, c(effect,rep(0,(ncol(infants)-length(acctrainstr)-1)),acctrainstr))
      
    }
}



###################### Leave one out Null #######################################################

# Calculate mean accuracy of prediction on final timepoint, via Leave One Out method of splitting Train-Test sets
LOO_multilogreg_Null <- function (infants,timepoint,end_timepoint,meta_file){
  
  logreg <- as.data.frame(matrix(NA,nrow = nrow(infants),ncol = 2))
  colnames(logreg)= c(paste('Cluster_at',end_timepoint,sep = '_'), timepoint)
  rownames(logreg)= rownames(infants)
  
  for (deigma in (rownames(infants))){
    
    logreg[deigma, 2] = ifelse(test = length(infants[deigma,timepoint]), yes = infants[deigma,timepoint], no = NA)
    logreg[deigma, 1] = ifelse(test = length(infants[deigma,as.character(end_timepoint)]), yes = infants[deigma,as.character(end_timepoint)], no = NA)
  
  }
  
  logreg = logreg[complete.cases(logreg),]
  
  acc <- c()
  tra <- c()
  form = as.formula(paste(paste('Cluster_at',end_timepoint,sep = '_'),'~.',sep = ''))
  for (j in 1:nrow(logreg)){
    trainset = logreg[1:nrow(logreg)!=j,]
    testset = logreg[j,]
    
    apotelesmata <- multinom(formula = form , data = trainset ,censored = T, model = T)
    
    giatotrain <- apotelesmata %>% predict(trainset)
    provlepseis <- apotelesmata %>% predict(testset)
    acc[j] = (sum(provlepseis == testset[,1])/ nrow(testset)) *100
    tra[j] = (sum(provlepseis == trainset[,1])/nrow(trainset))*100
    
  }
  
  return (c(mean(tra),mean(acc)))
}

###################### Stratified Train/Test splits Null ########################################

# Calculate accuracy of prediction on final timepoint, as a mean of 100 accuracies with different
# Train-Test stratified splits. Stratfication is performed on all metadata categories 
multilogreg_stratified_Null <- function(infants,pososto,timepoint,end_timepoint,meta_file){
  
  logreg <- as.data.frame(matrix(NA,nrow = nrow(infants),ncol = 2))
  colnames(logreg)= c(paste('Cluster_at',end_timepoint,sep = '_'),timepoint)
  rownames(logreg)= rownames(infants)
  
  for (deigma in (rownames(infants))){
  
    logreg[deigma, 2] = ifelse(test = length(infants[deigma,as.character(timepoint)]), yes = infants[deigma,as.character(timepoint)], no = NA)
    logreg[deigma, 1] = ifelse(test = length(infants[deigma,as.character(end_timepoint)]), yes = infants[deigma,as.character(end_timepoint)], no = NA)
    
  }
  logreg = logreg[complete.cases(logreg),]
  
  train_index <- createDataPartition(y = logreg[,paste('Cluster_at',end_timepoint,sep = '_')], p = pososto, list = F, times = 100, groups = max(apply(X = logreg, MARGIN = 2, FUN = function(x){length(unique(x))})) )
  
  acc <- c()
  tra <- c()
  paragontes <- list()
  fo = as.formula(paste(paste('Cluster_at',end_timepoint,sep = '_'),'~.',sep = ''))
  for (j in 1:ncol(train_index)){
    trainset = logreg[train_index[,j],]
    testset = logreg[-train_index[,j],]
    if (all(apply(X = logreg, MARGIN = 2, FUN = function(x){length(unique(x))})== apply(X = trainset, MARGIN = 2, FUN = function(x){length(unique(x))}))){
      
      apotelesmata <- multinom(formula = fo, data = trainset ,censored = T, model = T)
      
      giatotrain <- apotelesmata %>% predict(trainset)
      provlepseis <- apotelesmata %>% predict(testset)
      acc[j] = (sum(provlepseis == testset[,1])/ nrow(testset)) *100
      tra[j] = (sum(giatotrain == trainset[,1])/nrow(trainset))*100
    }
  }
  
  acc = na.omit(acc)
  tra = na.omit(tra)
  
  accteststr <- mean(acc)
  acctrainstr = mean(tra)
  trainsd <- sd(tra)
  testsd <- sd(acc)
  
  return (c(accteststr,acctrainstr, trainsd,testsd))
}

###################### Calculate optimal split percentage ############################################
# Calculate the optimal split percentage between .65 and .90 for predicting the final timepoint on the
# exact previous one to be used in all splits for all timepoints. The criterion is the best accuracy  
# as derived by the previous function with argument the different percentages.
bestpososto <- function (infants){
  olataaccuracies <- c()
  sepoiopososto <- c()
  for (poso in 65:90){
    poso = poso/100
    telika <- multilogreg_stratified_Null(infants = infants, pososto = poso, timepoint = tail(timepoints_to_perform,1), end_timepoint = end_timepoint,meta_file = files)
    accteststr = telika[1]
    acctrainstr = telika[2]
    trainsd = telika[3]
    testsd = telika[4]
    olataaccuracies= c(olataaccuracies,accteststr)
    sepoiopososto = c(sepoiopososto,poso)
  }
  
  return (sepoiopososto[which.max(olataaccuracies)])
}

# Assign the best percentage on a variable to use on the final calculations
to_kalitero_pososto <- bestpososto(infants = infants)

###################### Initializing variables Null ########################################
acctest <- c()
acctrain <-c()
accteststr <- c()
acctrainstr <- c()
trainsd <- c()
testsd <- c()
paragontes <- list()

###################### Calculate Accuracies on timepoints Null ############################

# Loop to create a model from every timepoint to the final timepoint as selected from the first loop
for (timepoint in timepoints_to_perform){
  
  LOOTeliko <- LOO_multilogreg_Null(infants = infants, timepoint = timepoint,end_timepoint = end_timepoint,meta_file = files)
  
  acctrain = c(acctrain,LOOTeliko[1])
  acctest = c(acctest,LOOTeliko[2])
  
  Splitteliko <- multilogreg_stratified_Null(infants = infants, pososto = to_kalitero_pososto,timepoint =  timepoint,meta_file =  files, end_timepoint = end_timepoint)
  
  
  accteststr <- c(accteststr,Splitteliko[1])
  acctrainstr <- c(acctrainstr, Splitteliko[2])
  trainsd <- c(trainsd, Splitteliko[3])
  testsd <- c(testsd, Splitteliko[4])
}  


###################### Save accuracies on the lists ######################################

# Adding to the lists in order to export the correct matrix
trainLOO = rbind (trainLOO, c('Null',rep(0,(ncol(infants)-length(acctrain)-1)),acctrain))
testLOO  = rbind (testLOO,  c('Null',rep(0,(ncol(infants)-length(acctest) -1)), acctest))
testStr  = rbind (testStr,  c('Null',rep(0,(ncol(infants)-length(accteststr) -1)), accteststr))
trainStr = rbind (trainStr, c('Null',rep(0,(ncol(infants)-length(acctrainstr)-1)), acctrainstr))






###################### Leave one out All #######################################
# Calculate mean accuracy of prediction on final timepoint, via Leave One Out method of splitting Train-Test sets
LOO_multilogreg <- function (infants,timepoint, end_timepoint){
  
  
  # Create the matrix with columns the state on the timepoint, the metadata and as a response the state on the last timepoint
  logreg <- as.data.frame(matrix(NA,nrow = nrow(infants),ncol = ncol(meta_file)))
  colnames(logreg)= c(paste('Cluster_at',end_timepoint,sep = '_'), timepoint , colnames(meta_file)[3:ncol(meta_file)])
  rownames(logreg)= rownames(infants)
  
  for (deigma in (rownames(infants))){
    subTab = meta_file[meta_file[,1] %in% deigma & meta_file[,2] == as.numeric(end_timepoint), 3:ncol(meta_file)]
    if (nrow(subTab) == 1){
      logreg[deigma, colnames(meta_file)[3:ncol(meta_file)]] = subTab[1,]
    }
    logreg[deigma, 2] = ifelse(test = length(infants[deigma,timepoint]), yes = infants[deigma,timepoint], no = NA)
    logreg[deigma, 1] = ifelse(test = length(infants[deigma,as.character(end_timepoint)]), yes = infants[deigma,as.character(end_timepoint)], no = NA)
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
    apotelesmata <- multinom(formula = fo , data = trainset ,censored = T, model = T)
    
    # Predict for trainset
    giatotrain <- apotelesmata %>% predict(trainset[,2:ncol(trainset)])
    # Predict for testset
    provlepseis <- apotelesmata %>% predict(testset[,2:ncol(testset)])
    # Calculate accuracies
    acc[j] = (sum(provlepseis == testset[,1])/ nrow(testset)) *100
    tra[j] = (sum(provlepseis == trainset[,1])/nrow(trainset))*100
    
  }
  
  return (c(mean(tra),mean(acc)))
}

###################### Stratified Train/Test splits All ########################

# Calculate accuracy of prediction on final timepoint, as a mean of 100 accuracies with different
# Train-Test stratified splits. Stratfication is performed on all metadata categories 
multilogreg_stratified <- function(infants,pososto,timepoint,times, end_timepoint){
  
  # Create the matrix with columns the state on the timepoint, the metadata and as a response the state on the last timepoint
  logreg <- as.data.frame(matrix(NA,nrow = nrow(infants),ncol = ncol(meta_file)))
  colnames(logreg)= c(paste('Cluster_at',end_timepoint,sep = '_'),timepoint,colnames(meta_file)[3:ncol(meta_file)])
  rownames(logreg)= rownames(infants)
  
  for (deigma in (rownames(infants))){
    subTab = meta_file[meta_file[,1] %in% deigma & meta_file[,2] == as.numeric(end_timepoint), 3:ncol(meta_file)]
    if (nrow(subTab) == 1){
      logreg[deigma, colnames(meta_file)[3:ncol(meta_file)]] = subTab[1,]
    }
    logreg[deigma, 2] = ifelse(test = length(infants[deigma,timepoint]), yes = infants[deigma,timepoint], no = NA)
    logreg[deigma, 1] = ifelse(test = length(infants[deigma,as.character(end_timepoint)]), yes = infants[deigma,as.character(end_timepoint)], no = NA)
  }
  # Remove NAs
  logreg = logreg[complete.cases(logreg),]
  
  # Create different partitions to split train and test sets
  train_index <- createDataPartition(y = logreg[,paste('Cluster_at',end_timepoint,sep = '_')], p = pososto, list = F, times = times, groups = max(apply(X = logreg, MARGIN = 2, FUN = function(x){length(unique(x))})) )
  
  # Set the response variable
  fo = as.formula(paste(paste('Cluster_at',end_timepoint,sep = '_'),'~.',sep = ''))
  
  acc <- c()
  tra <- c()
  paragontes <- list()
  for (j in 1:ncol(train_index)){
    
    # Split train and test sets
    trainset = logreg[train_index[,j],]
    testset = logreg[-train_index[,j],]
    if (all(apply(X = logreg, MARGIN = 2, FUN = function(x){length(unique(x))})== apply(X = trainset, MARGIN = 2, FUN = function(x){length(unique(x))}))){
      
      # Create a model
      apotelesmata <- multinom(formula = fo , data = trainset ,censored = T, model = T)
      # Predict for the trainset
      giatotrain <- apotelesmata %>% predict(trainset[,2:ncol(trainset)])
      # Predict for the testset
      provlepseis <- apotelesmata %>% predict(testset[,2:ncol(testset)])
      # Calculate accuracies
      acc[j] = (sum(provlepseis == testset[,1])/ nrow(testset)) *100
      tra[j] = (sum(giatotrain == trainset[,1])/nrow(trainset))*100
    }
  }
  # Remove NAs
  acc = na.omit(acc)
  tra = na.omit(tra)
  
  accteststr <- mean(acc)
  acctrainstr = mean(tra)
  trainsd <- sd(tra)
  testsd <- sd(acc)
  
  return (c(accteststr,acctrainstr, trainsd,testsd))
}

###################### Calculate optimal split percentage ######################
# Calculate the optimal split percentage between .65 and .90 for predicting the final timepoint on the
# exact previous one to be used in all splits for all timepoints. The criterion is the best accuracy  
# as derived by the previous function with argument the different percentages.
bestpososto <- function (infants){
  olataaccuracies <- c()
  sepoiopososto <- c()
  for (poso in 65:90){
    poso = poso/100
    telika <- multilogreg_stratified(end_timepoint = end_timepoint,infants = infants, pososto = poso, timepoint = rev(colnames(infants))[2], times = 100)
    accteststr = telika[1]
    acctrainstr = telika[2]
    trainsd = telika[3]
    testsd = telika[4]
    olataaccuracies= c(olataaccuracies,accteststr)
    sepoiopososto = c(sepoiopososto,poso)
  }
  
  return (sepoiopososto[which.max(olataaccuracies)])
}

# Assign the best percentage on a variable to use on the final calculations
to_kalitero_pososto <- bestpososto(infants = infants)


###################### Initializing variables All #############################

acctest <- c()
acctrain <-c()
accteststr <- c()
acctrainstr <- c()
trainsd <- c()
testsd <- c()
paragontes <- list()


###################### Calculate Accuracies on timepoints All #################

for (timepoint in timepoints_to_perform){
  
  
  LOOTeliko <- LOO_multilogreg(infants = infants, timepoint = timepoint,end_timepoint = end_timepoint)
  
  acctrain = c(acctrain,LOOTeliko[1])
  acctest = c(acctest,LOOTeliko[2])
  
  Splitteliko <- multilogreg_stratified(infants, to_kalitero_pososto, timepoint, times = 100,end_timepoint = end_timepoint)
  
  
  accteststr <- c(accteststr,Splitteliko[1])
  acctrainstr <- c(acctrainstr, Splitteliko[2])
  trainsd <- c(trainsd, Splitteliko[3])
  testsd <- c(testsd, Splitteliko[4])
}  

###################### Save accuracies on the lists ######################################
trainLOO = rbind (trainLOO,  c('All',rep(0,(ncol(infants)-length(acctrain)-1)),acctrain))
testLOO  = rbind (testLOO,   c('All',rep(0,(ncol(infants)-length(acctest)-1)),acctest))
testStr  = rbind (testStr,   c('All',rep(0,(ncol(infants)-length(accteststr)-1)),accteststr))
trainStr = rbind (trainStr,  c('All',rep(0,(ncol(infants)-length(acctrainstr)-1)),acctrainstr))

###################### Barplots of Accuracies ############################################
jpeg(filename = paste(dir_with_plots,paste('Accuracies of Models on Timepoint', end_timepoint,sep = ' '),sep = '/'))
barplot(rbind(acctest,accteststr, acctrain,acctrainstr), beside = T, ylim = c(40,100), xpd = F ,names = timepoints_to_perform, col= c('paleturquoise4','paleturquoise3', 'darkolivegreen2', 'darkolivegreen3'))
legend("topright", c("TestLOO","TestStratified", 'TrainLOO','TrainStratified'), fill=c('paleturquoise4','paleturquoise3', 'darkolivegreen2', 'darkolivegreen3'), cex = 1)
abline(h = 40.3, col = "white", lwd = 2, lty = 2)

dev.off()


}

###################### Write files of Accuracies ############################################
rownames(trainStr) = trainStr[,1]
rownames(trainLOO)= trainLOO[,1]
rownames(testLOO) = testLOO[,1]
rownames(testStr)= testStr[,1]
trainStr = trainStr[,2:ncol(trainStr)]
testStr = testStr[,2:ncol(testStr)]
trainLOO =trainLOO[,2:ncol(trainLOO)]
testLOO = testLOO[,2:ncol(testLOO)]

trainLOO = rbind(paste('From timepoint ',rev(colnames(infants)[1:ncol(infants)-1]),sep = ''),trainLOO)
testLOO  = rbind(paste('From timepoint ',rev(colnames(infants)[1:ncol(infants)-1]),sep = ''),testLOO)
trainStr = rbind(paste('From timepoint ',rev(colnames(infants)[1:ncol(infants)-1]),sep = ''),trainStr)
testStr  = rbind(paste('From timepoint ',rev(colnames(infants)[1:ncol(infants)-1]),sep = ''),testStr)
write.table(x = trainLOO, file = paste(dir_with_files,'Training Sets LOO.csv', sep = '/') , col.names = T, row.names = T)
write.table(x = testLOO,  file = paste(dir_with_files,'Test Sets LOO.csv',     sep = '/') , col.names = T, row.names = T)
write.table(x = trainStr, file = paste(dir_with_files,'Training Sets Splitted.csv', sep = '/') , col.names = T, row.names = T)
write.table(x = testStr,  file = paste(dir_with_files,'Test Sets Splitted.csv', sep = '/') , col.names = T, row.names = T)

write.csv(testLOO,file = paste(dir_with_files, 'TestLOO.csv',sep = '/'),row.names = T, col.names = c(paste('From timepoint ',rev(colnames(infants)[1:ncol(infants)-1]),sep = '')))
write.table(x = testLOO,file = paste(dir_with_files, 'TestLOO.csv',sep = '/'))
