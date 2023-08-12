get_total_effects <- function(weights,ix,end){
  
  TEs = rep(0,length(ix))
  for (i in seq_len(length(ix))){
    TEs[i] = walk_down(weights,ix[i])[end]
  }
  return(TEs)
}

walk_down <- function(weights,start){
  
  p = ncol(weights)
  total_effect = rep(0,p)
  beta = weights[start,] # take first step
  for (j in 1:p){
    total_effect = total_effect + beta # record into total effects
    beta = beta %*% weights # take next step
  }
  
  return(total_effect)
  
}

walk_down_recur <- function(weights, start, end, total_effect = 0){
  
  if (start==end) return(total_effect)
  
  ch = which(weights[start,] != 0)
  for (c in ch){
    total_effect = total_effect + walk_down_recur(weights,start,end, total_effect*weights[start,c]) # add the child and go through the child
  }
  
  return(total_effect)
  
}