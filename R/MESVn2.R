MESVn2 <- function(data,graph,Yi){
  
  p = ncol(data)-1
  S = matrix(0,nrow(data),p)
  E = matrix(0,nrow(data),p+1)
  
  # learn graph parameters
  direct_effects = matrix(0,p+1,p+1)
  for (i in 1:(p+1)){
    pa = which(graph[,i]!=0)
    if (length(pa)>0){
      mod = lm.fit(data[,pa,drop=FALSE],data[,i])
      direct_effects[pa,i] = mod$coefficients
      E[,i] = mod$residuals
    } else{
      E[,i]=data[,i]
    }
  }

  for (i in 1:p){
    
    S[,i] = get_total_effects(direct_effects,i,Yi)*E[,i]
    
  }
  
  return(S)
  
}

