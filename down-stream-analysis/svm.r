# train an svm classifier for every timepoint
svm.generate.classifiers <- function(otus, clusterings, percentage=0.8, kernel='linear'){

  # execute analysis per timepoint
  svms <- lapply(clusterings, function(clustering){
    
    # filter otus with timepoint
    otus.of.timepoint <- otus.all.records.for.timepoints(otus, clustering$timepoint)
    
    # normalize otus
    otus.normalized <- otus.normalize(otus.of.timepoint, 0)
    
    # train a classifier for every timepoint  
    svm.classifier <- svm.train.classifier(otus.normalized, clustering$best.clustering_by_ch, percentage, kernel)  
    
    result <- list('timepoint' = clustering$timepoint, 'svm' = svm.classifier)
    
    return(result)
  })
  
  return(svms)
}


svm.train.classifier <- function(data, classification, percentage, kernel){
  
  # if a timepoing has (n) clusters this list contains (n) classfication values 
  unique.classification.values <- unique(classification)
  
  # initiate training and validation data sets
  training.data.set <- array(dim=c(0,ncol(data)))
  training.classification.set <- array(dim=c(0,ncol(data)))
  validation.data.set <- array(dim=c(0,ncol(data)))
  validation.classification.set <- array(dim=c(0,ncol(data)))
  
  training.stats <- ''
  training.stats <- paste(training.stats, nrow(data), sep = '')
  training.stats <- paste(training.stats, ':', sep = '')
  
  # for every classification value create separate training and validation sets
  # Merge them all togenter and then train the SVMs
  for(classification.value in unique.classification.values){
  
    print(paste("Classification Value:", classification.value))
    
    classification.value.index <- as.vector(which(classification == classification.value))
    
    nrOfTrainingSamples <- round(percentage * length(classification.value.index), 0L)

    print(paste("Number of training samples:", nrOfTrainingSamples))

    training.stats <- paste(training.stats, '[', sep = '')
    training.stats <- paste(training.stats, classification.value, sep = '')
    training.stats <- paste(training.stats, ':', sep = '')
    
    if (nrOfTrainingSamples == 1){
      training.indeces <- classification.value.index
      validation.indeces <- classification.value.index
    }
    else{
      training.indeces <- sample(classification.value.index, nrOfTrainingSamples, replace = F)
      validation.indeces <- setdiff(classification.value.index, training.indeces)  
    }
    
    print("Training Indeces:")
    print(paste(training.indeces, sep = ' ,'))
    
    training.stats <- paste(training.stats, toString(length(training.indeces)), sep = '')
    training.stats <- paste(training.stats, '-', sep = '')
    training.stats <- paste(training.stats, toString(length(validation.indeces)), sep = '')

    training.data.set.for.specific.value <- data[training.indeces,]
    training.classification.set.for.specific.value <- classification[training.indeces]
   
    validation.data.set.for.specific.value <- data[validation.indeces,]
    validation.classification.set.for.specific.value <- classification[validation.indeces]   

    training.data.set <- rbind(training.data.set, training.data.set.for.specific.value)
    training.classification.set <- c(training.classification.set, training.classification.set.for.specific.value)
    
    validation.data.set <- rbind(validation.data.set, validation.data.set.for.specific.value)
    validation.classification.set <- c(validation.classification.set, validation.classification.set.for.specific.value)  
    
    training.stats <- paste(training.stats, ']', sep = '')
  }

  # train SVM with training set
  classifier <- svm(x = training.data.set, 
                    y = training.classification.set, 
                    type = "C-classification",
                    kernel = kernel, 
                    cost = 1, 
                    scale = F)
  
  # validate training
  prediction.on.validation.set <- svm.predict(classifier, validation.data.set)
  
  confusion.matrix <- table(prediction.on.validation.set, validation.classification.set)
  
  result <- list('classifier' = classifier,
                 'training.stats' = training.stats,
                 'confusion.matrix'= confusion.matrix)
  
  return(result)
}

svm.predict <- function(classifier, data){
  prediction <- predict(classifier, data)
  return(prediction)
}