setwd('Working_Chronos/Chronos_almost/')

source('Markov.R')
library (nnet)

ptm <- proc.time()

###############################################################################################################
############################ LEAVE ONE OUT ####################################################################
###############################################################################################################

# Calculate mean accuracy of prediction on final timepoint, via Leave One Out method of splitting Train-Test sets
LOO_multilogreg <- function (infants,timepoint){

  logreg <- as.data.frame(matrix(NA,nrow = nrow(infants),ncol = ncol(meta_file)))
  colnames(logreg)= c(paste('Cluster_at',tail(colnames(infants),1),sep = '_'), timepoint , colnames(meta_file)[3:ncol(meta_file)])
  rownames(logreg)= rownames(infants)
  
  for (deigma in (rownames(infants))){
    subTab = meta_file[meta_file[,1] %in% deigma & meta_file[,2] == 12, 3:ncol(meta_file)]
    if (nrow(subTab) == 1){
      logreg[deigma, colnames(meta_file)[3:ncol(meta_file)]] = subTab[1,]
    }
    logreg[deigma, 2] = ifelse(test = length(infants[deigma,tail(colnames(infants),2)[1]]), yes = infants[deigma,timepoint], no = NA)
    logreg[deigma, 1] = ifelse(test = length(infants[deigma,tail(colnames(infants),1)]), yes = infants[deigma,tail(colnames(infants),1)], no = NA)
  }
  
  logreg = logreg[complete.cases(logreg),]
  
  acc <- c()
  tra <- c()
  for (j in 1:nrow(logreg)){
    trainset = logreg[1:nrow(logreg)!=j,]
    testset = logreg[j,]
    
    apotelesmata <- multinom(formula = Cluster_at_24~. , data = trainset ,censored = T, model = T)
    
    giatotrain <- apotelesmata %>% predict(trainset[,2:ncol(trainset)])
    provlepseis <- apotelesmata %>% predict(testset[,2:ncol(testset)])
    acc[j] = (sum(provlepseis == testset[,1])/ nrow(testset)) *100
    tra[j] = (sum(provlepseis == trainset[,1])/nrow(trainset))*100
  
  }

 return (c(mean(tra),mean(acc)))
}


###############################################################################################################
############################ Stratified Train/Test splits #####################################################
###############################################################################################################

# Calculate accuracy of prediction on final timepoint, as a mean of 100 accuracies with different
# Train-Test stratified splits. Stratfication is performed on all metadata categories 
multilogreg_stratified <- function(infants,pososto,timepoint){
  
  logreg <- as.data.frame(matrix(NA,nrow = nrow(infants),ncol = ncol(meta_file)))
  colnames(logreg)= c(paste('Cluster_at',tail(colnames(infants),1),sep = '_'),timepoint,colnames(meta_file)[3:ncol(meta_file)])
  rownames(logreg)= rownames(infants)
  
  for (deigma in (rownames(infants))){
    subTab = meta_file[meta_file[,1] %in% deigma & meta_file[,2] == as.numeric(timepoint), 3:ncol(meta_file)]
    if (nrow(subTab) == 1){
      logreg[deigma, colnames(meta_file)[3:ncol(meta_file)]] = subTab[1,]
    }
    logreg[deigma, 2] = ifelse(test = length(infants[deigma,tail(colnames(infants),2)[1]]), yes = infants[deigma,timepoint], no = NA)
    logreg[deigma, 1] = ifelse(test = length(infants[deigma,tail(colnames(infants),1)]), yes = infants[deigma,tail(colnames(infants),1)], no = NA)
  }
  
  logreg = logreg[complete.cases(logreg),]

  train_index <- createDataPartition(y = logreg[,paste('Cluster_at',tail(colnames(infants),1),sep = '_')], p = pososto, list = F, times = 100, groups = max(apply(X = logreg, MARGIN = 2, FUN = function(x){length(unique(x))})) )
  
  acc <- c()
  tra <- c()
  paragontes <- list()
  for (j in 1:ncol(train_index)){
    trainset = logreg[train_index[,j],]
    testset = logreg[-train_index[,j],]
    if (all(apply(X = logreg, MARGIN = 2, FUN = function(x){length(unique(x))})== apply(X = trainset, MARGIN = 2, FUN = function(x){length(unique(x))}))){
      
      apotelesmata <- multinom(formula = Cluster_at_24 ~. , data = trainset ,censored = T, model = T)

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

# Calculate the optimal split percentage between .65 and .90 for predicting the final timepoint on the
# exact previous one to be used in all splits for all timepoints. The criterion is the best accuracy  
# as derived by the previous function with argument the different percentages.
bestpososto <- function (infants){
  olataaccuracies <- c()
  sepoiopososto <- c()
  for (poso in 65:90){
    poso = poso/100
    telika <- multilogreg_stratified(infants = infants, pososto = poso, timepoint = rev(colnames(infants))[2])
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
to_kalitero_pososto

###############################################################################################################
################################## Initializing variables #####################################################
###############################################################################################################

acctest <- c()
acctrain <-c()
accteststr <- c()
acctrainstr <- c()
trainsd <- c()
testsd <- c()
paragontes <- list()
###############################################################################################################
############################ Calculate Accuracies on timepoints  ##############################################
###############################################################################################################


for (timepoint in rev(colnames(infants))[2:ncol(infants)]){
  
  LOOTeliko <- LOO_multilogreg(infants = infants, timepoint = timepoint)
  
  acctrain = c(acctrain,LOOTeliko[1])
  acctest = c(acctest,LOOTeliko[2])
  
  Splitteliko <- multilogreg_stratified(infants, to_kalitero_pososto, timepoint)


  accteststr <- c(accteststr,Splitteliko[1])
  acctrainstr <- c(acctrainstr, Splitteliko[2])
  trainsd <- c(trainsd, Splitteliko[3])
  testsd <- c(testsd, Splitteliko[4])
}  

print (paste('Test set mean accuracy', round(mean(accteststr),2), sep = ' = '))
print (paste('Train set mean accuracy',round(mean(acctrainstr),2), sep = ' = '))

print (paste('Test set  LOO mean accuracy',round(mean(acctest),2), sep = ' = '))
print (paste('Train set LOO mean accuracy',round(mean(acctrain),2), sep = ' = '))

###############################################################################################################
####################################### Plotting  #############################################################
###############################################################################################################

#boxplot(acctest, acctrain ,names = c('TestSet Accuracy', 'TrainSet accuracy') , main = 'Accuracies Splits' ,outline = T, ylim = c(min(c(acctest,acctrain,accteststr,acctrainstr))-2,max(c(acctest,acctrain,accteststr,acctrainstr))+2))

#boxplot(accteststr, acctrainstr ,names = c('TestSet Accuracy', 'TrainSet accuracy') , main = 'Accuracies Splits' ,outline = T, ylim = c(min(c(accteststr,acctrainstr,acctest,acctrain))-2,max(c(accteststr,acctrainstr,acctest,acctrain))+2))

#boxplot(acctest,accteststr, acctrain, acctrainstr, names = c('TestSet LOO', 'TestSet Split', 'TrainSet LOO', 'TrainSet Split'), outline = T, ylim = c(min(c(accteststr,acctrainstr,acctest,acctrain))-2,max(c(accteststr,acctrainstr,acctest,acctrain))+2))


xrwma <- c('paleturquoise4','paleturquoise3', 'mediumturquoise', 'darkturquoise')
for (i in 1:(ncol(infants)-1)){
  jpeg(filename =paste(dir_with_plots, paste(paste('Multinomial Logistic Regression Accuracies',colnames(infants[i]),sep = ' of Timepoint '),'jpeg',sep = ' .'),sep = '/'))
  barplot(cbind(acctest,accteststr, acctrain, acctrainstr)[i,] , name = c('Test LOO','Test Splits', 'Train LOO','Train Split') ,col= xrwma, main =paste('Mean Accuracy on Timepoint',rev(colnames(infants)[2:ncol(infants)])[i], sep = ' ') ,ylim = c(0,100) )
  dev.off()
}

proc.time() - ptm

###############################################################################################################
####################################### TO DO #################################################################
###############################################################################################################

# Multinomial Logistic Regression on INFANTS (to check if a profiles on multiple timepoints can be predictive)

# Find how to use formula globally


acctest
acctrain
accteststr
acctrainstr
