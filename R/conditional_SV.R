conditional_SV <- function(X,Y,nsamp=1000){
  require(earth)
  
  S = X
  p = ncol(X)
  for (i in 1:p ){
    # ix = setdiff(1:p,i)
    for (j in 1:nsamp){
      perm = sample(1:p,p,replace=FALSE)
      k = which(perm==i)
      if (k>1){
        W = perm[1:(k-1)]
      } else{
        W = NULL
      }
      
      psi1 = earth(x=X[,c(W,i),drop=FALSE],y=Y)$fitted.values
      
      if (k == 1){
        psi2 = mean(Y)
      } else{
        psi2 = earth(x=X[,W,drop=FALSE],y=Y)$fitted.values
        
      }
      
      gamma = psi1 - psi2
      
      S[,i] = S[,i] + gamma
      
      
    }
    
    S[,i] = S[,i]/nsamp
    
  }
  
  return(S)
}