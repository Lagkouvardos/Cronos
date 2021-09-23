#setwd('Working_Chronos/Chronos_almost/')

source('Markov.R')
library (nnet)


akjghfdsa <- formula(x = colnames(logreg)[1])
logreg
###############################################################################################################
############################ LEAVE ONE OUT ####################################################################
###############################################################################################################

acctrain <- c()
acctest <- c()

for (i in 1:ncol(infants)-1){
  logreg <- as.data.frame(matrix(NA,nrow = nrow(infants),ncol = ncol(meta_file)))
  colnames(logreg)= c(paste('Cluster_at',tail(colnames(infants),1),sep = '_'), colnames(infants)[ncol(infants)-i] , colnames(meta_file)[3:ncol(meta_file)])
  rownames(logreg)= rownames(infants)
  
  for (deigma in (rownames(infants))){
    subTab = meta_file[meta_file[,1] %in% deigma & meta_file[,2] == 12, 3:ncol(meta_file)]
    if (nrow(subTab) == 1){
      logreg[deigma, colnames(meta_file)[3:ncol(meta_file)]] = subTab[1,]
    }
    logreg[deigma, 2] = ifelse(test = length(infants[deigma,tail(colnames(infants),2)[1]]), yes = infants[deigma,colnames(infants)[ncol(infants)-i]], no = NA)
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
  
  acctest[i] <- mean(acc)
  acctrain[i] = mean(tra)

}


###############################################################################################################
############################ Stratified Train/Test splits #####################################################
###############################################################################################################


accteststr <- c()
acctrainstr <- c()
trainsd <- c()
testsd <- c()
for (i in 1:ncol(infants)-1){
  logreg <- as.data.frame(matrix(NA,nrow = nrow(infants),ncol = ncol(meta_file)))
  colnames(logreg)= c(paste('Cluster_at',tail(colnames(infants),1),sep = '_'),colnames(infants)[ncol(infants)-i],colnames(meta_file)[3:ncol(meta_file)])
  rownames(logreg)= rownames(infants)
  
  for (deigma in (rownames(infants))){
    subTab = meta_file[meta_file[,1] %in% deigma & meta_file[,2] == 12, 3:ncol(meta_file)]
    if (nrow(subTab) == 1){
      logreg[deigma, colnames(meta_file)[3:ncol(meta_file)]] = subTab[1,]
    }
    logreg[deigma, 2] = ifelse(test = length(infants[deigma,tail(colnames(infants),2)[1]]), yes = infants[deigma,colnames(infants)[ncol(infants)-i]], no = NA)
    logreg[deigma, 1] = ifelse(test = length(infants[deigma,tail(colnames(infants),1)]), yes = infants[deigma,tail(colnames(infants),1)], no = NA)
  }
  
  logreg = logreg[complete.cases(logreg),]

  train_index <- createDataPartition(y = logreg[,paste('Cluster_at',tail(colnames(infants),1),sep = '_')], p = 0.85, list = F, times = 100, groups = max(apply(X = logreg, MARGIN = 2, FUN = function(x){length(unique(x))})) )
  
  acc <- c()
  tra <- c()
  
  for (j in 1:ncol(train_index)){
    trainset = logreg[train_index[,j],]
    testset = logreg[-train_index[,j],]
    if (all(apply(X = logreg, MARGIN = 2, FUN = function(x){length(unique(x))})== apply(X = trainset, MARGIN = 2, FUN = function(x){length(unique(x))}))){
      apotelesmata <- multinom(formula = Cluster_at_24 ~. , data = trainset ,censored = T, model = T)
      #paragontes[[i]] <- kati_kalo$coefficients
      giatotrain <- apotelesmata %>% predict(trainset[,2:ncol(trainset)])
      provlepseis <- apotelesmata %>% predict(testset[,2:ncol(testset)])
      acc[j] = (sum(provlepseis == testset[,1])/ nrow(testset)) *100
      tra[j] = (sum(giatotrain == trainset[,1])/nrow(trainset))*100
    }
    else {
      next
    }
  }
  
  acc = na.omit(acc)
  tra = na.omit(tra)
  
  accteststr[i] <- mean(acc)
  acctrainstr[i] = mean(tra)
  trainsd[i] <- sd(tra)
  testsd[i] <- sd(acc)

}

print (paste('Test set mean accuracy', round(mean(accteststr),2), sep = ' = '))
print (paste('Train set mean accuracy',round(mean(acctrainstr),2), sep = ' = '))

print (paste('Test set  LOO mean accuracy',mean(acctest), sep = ' = '))
print (paste('Train set LOO mean accuracy',mean(acctrain), sep = ' = '))

###############################################################################################################
####################################### Plotting  #############################################################
###############################################################################################################

#boxplot(acctest, acctrain ,names = c('TestSet Accuracy', 'TrainSet accuracy') , main = 'Accuracies Splits' ,outline = T, ylim = c(min(c(acctest,acctrain,accteststr,acctrainstr))-2,max(c(acctest,acctrain,accteststr,acctrainstr))+2))

#boxplot(accteststr, acctrainstr ,names = c('TestSet Accuracy', 'TrainSet accuracy') , main = 'Accuracies Splits' ,outline = T, ylim = c(min(c(accteststr,acctrainstr,acctest,acctrain))-2,max(c(accteststr,acctrainstr,acctest,acctrain))+2))


boxplot(acctest,accteststr, acctrain, acctrainstr, names = c('TestSet LOO', 'TestSet Split', 'TrainSet LOO', 'TrainSet Split'), outline = T, ylim = c(min(c(accteststr,acctrainstr,acctest,acctrain))-2,max(c(accteststr,acctrainstr,acctest,acctrain))+2))

###############################################################################################################
####################################### TO DO #################################################################
###############################################################################################################


