do_SVn2 <- function(data,graph,Yi,exc=NULL,nsamp=100){
  require(earth)
  
  X = data[,-Yi]
  Y = data[,Yi]
  
  p = ncol(X)
  S = matrix(0,nrow(X),p)

  direct_effects = matrix(0,p+1,p+1)
  for (i in 1:(p+1)){
    pa = which(graph[,i]!=0)
    if (length(pa)>0){
      direct_effects[pa,i] = lm.fit(data[,pa,drop=FALSE],data[,i])$coefficients
    } 
  }
  
  poss = setdiff(1:p,exc)
  for (i in poss ){
    for (j in 1:nsamp){
      if (length(poss>1)){
        perm = sample(poss,length(poss),replace=FALSE) 
      } else{
        perm = i
      }
      k = which(perm==i)
      if (k>1){
        W = perm[1:(k-1)]
      } else{
        W = NULL
      }
      
      direct_effectsT = direct_effects
      direct_effectsT[,c(W,i)] = 0
      psi1 = data[,c(W,i),drop=FALSE] %*% as.matrix(get_total_effects(direct_effectsT,c(W,i),Yi))
      
      if (length(W)==0){
        psi2 = mean(Y)
      } else{
        direct_effectsT = direct_effects
        direct_effectsT[,W] = 0
        psi2 = data[,W,drop=FALSE] %*% as.matrix(get_total_effects(direct_effectsT,W,Yi))
      }
     
      S[,i] = S[,i] + (psi1 - psi2)
      
    }
    
    S[,i] = S[,i]/nsamp
    
  }
  
  return(S)
}
