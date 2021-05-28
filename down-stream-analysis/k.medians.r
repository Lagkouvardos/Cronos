source('base.r')

# calculate manhattan distance between vectors p, q
manhattan.distance <- function(p, q){
  return (sum(abs(p - q)))
}

# calculate median of input vector
kmedians.medians.of.columns <- function(X){
  
  medians = list()
  nr_of_columns <- ncol(X)

  for (i in 1:nr_of_columns){
    medians[[i]] = median(X[,i], na.rm = F)
  }
  
  return(unlist(medians))
}

# initialize centroids by randomly selecting vectors from the dataset
kmedians.init.centroids.with.random.selection <- function(X, nr_of_centroids){
  
  centroids <- list()
  nr_of_dimensions = ncol(X)
  
  random_indexes <- sample.int(nrow(X), nr_of_centroids)
  
  for (i in 1:nr_of_centroids){
    index = random_indexes[i]
    centroids[[i]] <- as.vector(X[index, ])
  }
  
  return (centroids)
}

# initialize centroids by calculating k-means clusters
kmedians.init.centroids.with.kmeans <- function(X, nr_of_centroids){

  centroids <- list()
  nr_of_dimensions = ncol(X)  
  kmeans_result <- kmeans(X, nr_of_centroids, nstart = 25)
  kmeans_centroids <- as.matrix(kmeans_result["centers"][[1]])
  
  for (i in 1:nrow(kmeans_centroids)){
    centroids[[i]] <- as.vector(kmeans_centroids[i, ])
  }
  
  return (centroids)
}

# initialize centroids
kmedians.init.centroids <- function(X, nr_of_centroids, type){
  
  if (type == 'k-means')
    return(kmedians.init.centroids.with.kmeans(X, nr_of_centroids))
    
  return(kmedians.init.centroids.with.random.selection(X, nr_of_centroids))
}

# calculate centroid that is nearest to data_vector
kmedians.nearest.centroid <- function(centroids, data_vector){
  
  nr_of_centroids = length(centroids)
  distance_from_centroids = list()
  
  for (k in 1:nr_of_centroids){
    centroid = centroids[[k]]
    distance_from_centroids[[k]] = manhattan.distance(data_vector, centroid)
  }
  
  min_distance = min(unlist(distance_from_centroids))
  centroid_index = which.min(unlist(distance_from_centroids))
  
  return (list(centroid_index, min_distance))
}

# allocate data vectors to centroids
kmedians.allocate.to.centroids <- function(X, centroids){
  
  nr_of_items = nrow(X)
  centroid_allocations = matrix(nrow = nr_of_items, ncol = 1)
  
  for (i in 1:nr_of_items){
    data_vector <- X[i,]
    nearest_centroid <- kmedians.nearest.centroid(centroids, data_vector)
    centroid_allocations[i,1] = nearest_centroid[[1]]
  }
  
  return(centroid_allocations)
}

# update centroids based on median of latest allocations
kmedians.update.centroids <- function(X, centroids, centroid_allocations){
  
  new_centroids <- list()
  nr_of_centroids <- length(centroids)
  
  for (k in 1:nr_of_centroids){
    item_indexes_allocated_to_centroid <- which(centroid_allocations[, 1] == k)
    items_allocated_to_centroid <- as.matrix(X[item_indexes_allocated_to_centroid, ])

    if (nrow(items_allocated_to_centroid) == 0){
      new_centroids[[k]] <- integer(ncol(X))
    }
    else{
      items_medians <- kmedians.medians.of.columns(items_allocated_to_centroid)
      new_centroids[[k]] <- items_medians
    }
  }
  
  return(new_centroids)
}

# perform kmedians clustering on input data
kmedians <- function(X, centers, max_iterations = 100){
  
  centroids <- kmedians.init.centroids(X, centers, 'k-means')

  for (i in 1:max_iterations){

    centroid_allocations = kmedians.allocate.to.centroids(X, centroids)
    new_centroids <- kmedians.update.centroids(X, centroids, centroid_allocations)
    print(paste('Update centroids - iteration :', i, sep = " "))

    no_change_in_centroids = vectors.are.equal(unlist(centroids), unlist(new_centroids))
    if (no_change_in_centroids) { break }

    centroids = new_centroids
  }
  
  rownames(centroid_allocations) <- row.names(X)
  return(centroid_allocations)
}