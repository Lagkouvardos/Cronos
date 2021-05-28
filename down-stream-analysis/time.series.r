time.series.cluster.timepoint <- function(otus.normalized, distances, timepoint){
  
  # cluster distances using k-medoids
  max_nr_of_clusters <- nrow(otus.normalized) - 1
  clustering.instances <- lapply(2:max_nr_of_clusters, function(nr_of_clusters){
    kmedoids <- pam(distances, nr_of_clusters, diss = T)
    calinski_harabasz_index <- calinhara(otus.normalized, kmedoids$clustering, nr_of_clusters)
    silhouette_index <- kmedoids$silinfo$avg.width
    list('clustering' = kmedoids$clustering, 'ch_benchmark' = calinski_harabasz_index, 'slh_benchmark' = silhouette_index)
  })
  
  ch_benchmark_list <- unlist(lapply(clustering.instances, function(instance){instance$ch_benchmark}))
  slh_benchmark_list <- unlist(lapply(clustering.instances, function(instance){instance$slh_benchmark}))
  
  best_clustering_index_by_slh <- which(slh_benchmark_list == max(slh_benchmark_list))
  maxL <- (length(ch_benchmark_list) / 2) + 1
  best_clustering_index_by_ch <- which(ch_benchmark_list == max(ch_benchmark_list[1:maxL]))
  
  best_nr_of_clusters_by_slh <- best_clustering_index_by_slh + 1
  best_nr_of_clusters_by_ch <- best_clustering_index_by_ch + 1
  
  best_clustering_by_slh <- clustering.instances[[best_clustering_index_by_slh]]$clustering
  best_clustering_by_ch <- clustering.instances[[best_clustering_index_by_ch]]$clustering

  # build result
  result <- list('best.clustering_by_slh' = best_clustering_by_slh, 
                 'best.clustering_by_ch' = best_clustering_by_ch,
                 'best.nr.of.clusters_by_slh' = best_nr_of_clusters_by_slh,
                 'best.nr.of.clusters_by_ch' = best_nr_of_clusters_by_ch,
                 'ch_benchmark_list' = ch_benchmark_list,
                 'slh_benchmark_list' = slh_benchmark_list,
                 'timepoint' = timepoint,
                 'distances' = distances)
  
  return (result)
}

time.series.generate.clusterings <- function(otus, otus.tree){
  
  # execute analysis per timepoint
  clusterings <- lapply(unique.timepoints, function(timepoint){
    
    # filter otus with timepoint
    otus.of.timepoints <- otus.all.records.for.timepoints(otus, timepoint)
    
    # normalize otus
    otus.normalized <- otus.normalize(otus.of.timepoints, 0)
    
    # calcuate distances
    distances <- distances.calculate(otus.normalized, otus.tree)
    
    # cluster data
    time.series.cluster.timepoint(otus.normalized, distances, timepoint)
  })
  
  return(clusterings)
}

time.series.generate <- function(clusterings, samples){
  
  individuals <- unlist(lapply(samples, function(i){substr(i, start = 2, stop = 4)}))
  timepoints <- unlist(lapply(samples, function(i){substr(i, start = 5, stop = 7)}))
  
  time.series <- matrix(nrow = length(unique(individuals)), ncol = length(unique(timepoints)))
  row.names(time.series) <- unique(individuals)
  colnames(time.series) <- unique(timepoints)
  
  for (clustering in clusterings)
  {
    for (individual in unique(individuals))
    {
      timepoint <- clustering$timepoint
      sample = paste('X', individual, timepoint, sep = '')
      time.series[individual, timepoint] = clustering$best.clustering_by_ch[sample]
    }
  }
  
  return(time.series)
}