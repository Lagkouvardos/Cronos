vectors.are.equal <- function(a, b){
  
  lenght_of_a <- length(a)
  lenght_of_b <- length(b)
  
  if (lenght_of_a != lenght_of_b) {return (FALSE)}
  
  for (i in 1:lenght_of_a){
    if(a[i] != b[i]) {return (FALSE)}
  }
  
  return (TRUE)
}

select.item.with <- function(X, name, value){
  
  items <- lapply(X, function(item){ if(item[name] == value) item else NA})
  result <- items[which(!is.na(items))]
  return(result)
}