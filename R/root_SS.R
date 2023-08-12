root_SS <- function(data,truth){

  data = as.matrix(data)
  
  cl <- hclust.vector(data, method="ward")
  
  SS = rep(0,10)
  # SS_t = rep(0,10)
  for (c in 1:10){
    ix = cutree(cl, k = c)
    for (cc in 1:c){
      for (d in 1:ncol(data))
        SS[c] = SS[c] + sum(( truth[which(ix==cc),d] - mean(data[which(ix==cc),d]) )^2)
        # SS_t[c] = SS_t[c] + sum(( truth[which(ix==cc),d] - mean(truth[which(ix==cc),d]) )^2)
    }
  }

  return(SS/nrow(data))

}