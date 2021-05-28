# plot timepoints
plots.timepoints <- function(x, y, samples, timepoints){
  fig_by_month <- plot_ly(x = x, y = y,
          type = "scatter",
          color = timepoints, colors = brewer.pal(n = length(timepoints), name = "Set1"),
          mode = "markers") %>% 
          layout(title = 'Timepoint')
  fig_by_month
}

# plot samples
plots.samples <- function(x, y, samples, timepoints){
	fig_by_sample <- plot_ly(x = x, y = y,
							 color = samples, colors = brewer.pal(n = length(samples), name = "Set1"),
							 text = timepoints,
							 type = "scatter", mode = "markers") %>% 
			layout(title = 'Samples') %>% 
			add_text(textposition = "top right")
	fig_by_sample
}

# plot classification
plots.classification <- function(x, y, centroid_allocations){
  fig_by_classification <- plot_ly(x = x, y = y,
           color = as.vector(centroid_allocations), colors = brewer.pal(n = length(centroid_allocations), name = "Set1"),
           type = "scatter", mode = "markers")
   fig_by_classification <- fig_by_classification %>% layout(title = 'Classification')
   fig_by_classification
}

#plot specific individual
plots.individual <- function(x, y, samples, timepoints, individual){
  
  individual.index <- which(startsWith(samples, paste("X", individual, sep = '')))
  
  x.individual <- x[individual.index]
  y.individual <- y[individual.index]
  
  individual.timepoints <- timepoints[individual.index]
  individual.samples <- samples[individual.index]
  
  fig_by_sample <- plot_ly(x = x.individual, y = y.individual,
                           text = individual.timepoints,
                           type = "scatter", mode = "markers") %>% 
    layout(title = individual) %>% 
    layout(showlegend = FALSE) %>% 
    add_text(textposition = "top right")
  fig_by_sample
}
  
plots.network <- function(transition.matrix, time.point.list){
  
  g <- graph_from_adjacency_matrix(transition.matrix, weighted = T, mode = 'directed')
  edge.list <- get.data.frame(g, what = 'edges')
  node.list <- get.data.frame(g, what = 'vertices')

  x <- c()
  y <- c()
  tm.rownames <- rownames(transition.matrix)
  for (rowname in tm.rownames){

    timepoint <- unlist(strsplit(rowname, ".", fixed = T))[1]
    cluster <- unlist(strsplit(rowname, ".", fixed = T))[2]
    x <- c(x, which(timepoint == time.point.list))
    
    rownames.of.timepoint <- tm.rownames[str_starts(tm.rownames, timepoint)]
    y <- c(y, which(endsWith(rownames.of.timepoint, cluster)))
  }
  
  node.list$x <- x
  node.list$y <- y
  
  new_g <- graph_from_data_frame(d = edge.list, directed = T, vertices = node.list)
  plot(new_g,
       vertex.size=30, 
       edge.label=round(edge.list$weight, 2), 
       edge.width=round(edge.list$weight * 5, 2),
       xlim=c(0,0))
}

plots.multi.timepoint.clustering <- function(mds, clusterings){

  rownames.of.x <- rownames(mds)
  
  # adjust clustering values, so that plot of clustering per timepoint in readable ... 
  merged.clustering.by.ch <- array(dim = nrow(mds))
  
  index = 1
  for(i in rownames.of.x){
    sample <- substr(i, start = 2, stop = 4)
    timepoint <- substr(i, start = 5, stop = 7)
    
      clustering <- NULL
      for(j in clusterings)
        if(j$timepoint == timepoint) clustering = j
      
      value <- as.integer(clustering$best.clustering_by_ch[i])
      merged.clustering.by.ch[index] = paste(timepoint, value, sep='-')
      index = index + 1
    
  }
  
  fig_by_classification <- plot_ly(x = mds[,1], y = mds[,2], type = "scatter", mode = "markers" 
                                   ,text = as.vector(merged.clustering.by.ch)
                                   ,color = as.vector(merged.clustering.by.ch)
                                   ,colors = brewer.pal(n = length(merged.clustering.by.ch), name = "Set1")
                                   )
  fig_by_classification <- fig_by_classification %>% add_text(textposition = "top right")
  fig_by_classification <- fig_by_classification %>% layout(title = 'Classification')
  fig_by_classification
}
  