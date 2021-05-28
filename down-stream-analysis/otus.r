otus.all.records.for.sample <- function(X, sample_id){
  row_names <- row.names(X)
  indeces <- which(startsWith(row_names, paste("X", sample_id, sep = '')))
  return(X[indeces, ])
}

otus.all.records.for.timepoints <- function(otus, timepoint){
  row_names <- row.names(otus)
  indeces <- which(substr(samples, start = 5, stop = 7) %in% timepoint)
  return(otus[indeces, ])
}

otus.normalize <- function(otus, otu_threashold){
  
  # calculate otu percentages per sample
  otus_sums_by_sample <- rowSums(otus)
  otus_percentages <- otus / otus_sums_by_sample
  
  # filter otus. if an otu is below threashold for all samples, exclude it
  otus_to_keep <- colSums(otus_percentages > otu_threashold) > 0
  otus_percentages_filtered <- otus_percentages[, otus_to_keep]
  
  # normalize otus
  otus_norm <- min(rowSums(otus)) * otus_percentages_filtered
  #mode(otus_norm) <- "integer"
  
  return(otus_norm)
}

otus.load <- function(data.folder, otus.filename, otus.tree.filename){
  
  # read otus file
  otus <- t(read.delim2(paste(data.folder, otus.filename, sep = ''), header=T, sep="\t", row.names = 1))
  
  # remove taxonomy row from otu table
  otus.contains.taxonomy.row <- tail(rownames(otus), n=1) == "taxonomy"
  if (otus.contains.taxonomy.row) {otus <- otus[1:nrow(otus)-1,]}
  
  # make all strings integers
  mode(otus) <- "integer"
  
  # order otus based on rownames (aka sample names)
  otus <- otus[ order(row.names(otus)), ]
  
  # read the tree
  otus.tree <- read.tree(paste(data.folder, otus.tree.filename, sep = ''))

  result <- list('otus' = otus, 'otus.tree' = otus.tree)
  return(result)
}