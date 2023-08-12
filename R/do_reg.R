do_reg <- function(data,graph,Yi,exc=NULL){
  
  p = ncol(data)-1
  S = matrix(0,nrow(data),p)
  
  # learn graph parameters
  direct_effects = matrix(0,p+1,p+1)
  for (i in 1:(p+1)){
    pa = which(graph[,i]!=0)
    if (length(pa)>0){
      direct_effects[pa,i] = lm.fit(data[,pa,drop=FALSE],data[,i])$coefficients
    } 
  }
  
  for (i in setdiff(1:p,exc)){
   
      S[,i] = get_total_effects(direct_effects,i,Yi)*data[,i]

  }
  
  return(S)
  
}

