markov.generate.state.series <- function(time.series){
  
  column.names <- colnames(time.series)
  
  state.series <- time.series
  for(i in 1:nrow(time.series)){
    for(j in 1:ncol(time.series)){
      if(is.na(time.series[i,j])){
        #TODO
      }
      else{
        state.series[i,j] <- paste(column.names[j], time.series[i,j], sep = '.')  
      }
    }  
  }
  
  return(state.series)
}

markov.transition.probability <- function(state.series, from.state, to.state){
  
  from.timepoint <- unlist(strsplit(from.state, ".", fixed = T))[1]
  to.timepoint <- unlist(strsplit(to.state, ".", fixed = T))[1]
  
  timepoint.states <- state.series[, from.timepoint]
  timepoint.states <- timepoint.states[!is.na(timepoint.states)]
  from.state.count <- count(timepoint.states == from.state)
  
  to.state.count <- 0
  from.state.count.reduction <- 0
  
  for (i in 1:nrow(state.series)){

    state.serie <- state.series[i, ]    
    from.serie.state <- state.serie[from.timepoint]
    to.serie.state <- state.serie[to.timepoint]
    
    to.is.empty <- is.null(to.serie.state) || is.na(to.serie.state)
    from.is.empty <- is.null(from.serie.state) || is.na(from.serie.state)
      
    is.empty <- to.is.empty || from.is.empty
      
    if (!is.empty && to.serie.state == to.state && from.serie.state == from.state){
      to.state.count <- to.state.count + 1  
    }
    
    if (to.is.empty && !from.is.empty && from.serie.state == from.state){
      from.state.count.reduction <- from.state.count.reduction + 1
    }
  }
  
  return(to.state.count / (from.state.count - from.state.count.reduction))
}

markov.transition.matrix.for.timepoint.pair <- function(state.series, from.timepoint.states, to.timepoint.states){
  
  all.pair.states <- c(from.timepoint.states, to.timepoint.states)
  pair.transition.matrix <- matrix(0L, nrow = length(all.pair.states), ncol = length(all.pair.states))
  rownames(pair.transition.matrix) <- all.pair.states
  colnames(pair.transition.matrix) <- all.pair.states

  for (from.state in from.timepoint.states){
    for (to.state in to.timepoint.states){
      i = which(all.pair.states == from.state)
      j = which(all.pair.states == to.state)
      pair.transition.matrix[i,j] <- markov.transition.probability(state.series, from.state, to.state)
    }  
  }
  
  return(pair.transition.matrix)
}

markov.merge.pair.transition.matrix <- function(transition.matrix, pair.transition.matrix){
  
  for(i in 1:nrow(pair.transition.matrix)){
    for(j in 1:ncol(pair.transition.matrix)){
      
      row.name <- rownames(pair.transition.matrix)[i] 
      column.name <- colnames(pair.transition.matrix)[j]
      
      i_new <- which(rownames(transition.matrix) == row.name)
      j_new <- which(colnames(transition.matrix) == column.name)
      
      transition.matrix[i_new, j_new] <- pair.transition.matrix[i, j]
    }  
  }
  
  return(transition.matrix)
}

markov.transition.matrix <- function(time.series, time.point.list = NULL){
  
  state.series <- markov.generate.state.series(time.series)
  states.per.timepoint <- lapply(1:ncol(state.series), function(i){unique(state.series[,i])})
  states.per.timepoint <- lapply(states.per.timepoint, function(i){i[!is.na(i)]})
  
  if (is.na(time.point.list) || is.null(time.point.list)){
    time.point.list <- colnames(time.series)
  }
  
  states.for.selected.timepoints <- c()
  for(time.point in time.point.list){
    time.point.index = which(colnames(time.series) == time.point)
    states.for.selected.timepoints <- c(states.for.selected.timepoints, states.per.timepoint[time.point.index])
  }

  transition.matrix <- matrix(0L, nrow = length(unlist(states.for.selected.timepoints)), ncol = length(unlist(states.for.selected.timepoints)))
  rownames(transition.matrix) <- unlist(states.for.selected.timepoints)
  colnames(transition.matrix) <- unlist(states.for.selected.timepoints)
  
  for (i in 1:length(states.for.selected.timepoints)){

    from.timepoint.states <- unlist(states.for.selected.timepoints[i])
    to.timepoint.states <- unlist(states.for.selected.timepoints[i+1])
    
    if (!is.null(to.timepoint.states)){
      pair.transition.matrix <- markov.transition.matrix.for.timepoint.pair(state.series, from.timepoint.states, to.timepoint.states)  
      transition.matrix <- markov.merge.pair.transition.matrix(transition.matrix, pair.transition.matrix)
    }
  }
  
  return(transition.matrix)
}