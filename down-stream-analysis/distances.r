distances.calculate <- function(otus.normalized, otus.tree, distance.type = "unifrac"){
 
  distance <- NULL
  if (distance.type == "unifrac"){
    unifracs <- GUniFrac(otus.normalized, midpoint(otus.tree), alpha = c(0.0, 0.5, 1.0))$unifracs
    distance <- as.dist(unifracs[, , "d_0.5"])
  }
  else if (distance.type == "euclidian") {
    distance <- dist(otus.normalized,  method = 'euclidian', upper = T)
  }
  else if (distance.type == "manhattan") {
    distance <- dist(otus.normalized,  method = 'manhattan', upper = T)
  } 
  
  if (is.null(distance) || is.na(distance)) {
    stop("time.series.calculate.distances() function called with invalid distance type")
    return(NA)
  }
  
  return(distance)
}