environment.list.of.libraries <- function(){
	return(c('ape', 
	    'cluster', 
	    'fpc', 
	    'phangorn', 
	    'GUniFrac', 
	    'plotly', 
	    'RColorBrewer', 
	    'stringr',
	    'igraph',
	    'e1071',
	    'readxl',
	    'randomForest',
	    'bnlearn'))
}

environment.load.sources <- function(){
  
  # load .r code from local files
  source('base.r')
  source('environment.r')
  source('otus.r')
  source('plots.r')
  source('time.series.r')
  source('markov.r')
  source('svm.r')
  source('distances.r')
  source('random.forest.r')
}

environment.install.missing.packages <- function(){
  
  required.packages <- as.vector(environment.list.of.libraries())

  already.installed.packages <- installed.packages()[,"Package"]
  
  matching.packages <- match(required.packages, already.installed.packages)
  
  missing.packages <- required.packages[which(is.na(matching.packages))]
  
  if (length(missing.packages) > 0) { install.packages(missing.packages, dependencies = T)}
}

environment.load.packages <- function(){
  
  libraryDeclaratioText <- " "
  for(entry in environment.list.of.libraries()){
    libraryDeclaratioText <- paste(libraryDeclaratioText, "library('", sep = "")
    libraryDeclaratioText <- paste(libraryDeclaratioText, entry, sep = "")
    libraryDeclaratioText <- paste(libraryDeclaratioText, "')", sep = "")
    libraryDeclaratioText <- paste(libraryDeclaratioText, "\n", sep = "")
  }

  eval(parse(text=libraryDeclaratioText))
}

environment.load.dependencies <- function(){
  
  environment.install.missing.packages()
  environment.load.packages()
  environment.load.sources()
}

environment.cleanup <- function(){
  # cleanup
  cat("\014")  
  rm(list=ls(all=TRUE))
}

environment.start <- function(){
  environment.load.dependencies()
  environment.cleanup()
}