marg_SVn2 <- function(data,graph,Yi,nsamp=100){
  #linear
  X = data[,-Yi]
  Y = data[,Yi]
  
  p = ncol(X)
  S = matrix(0,nrow(X),p)
  
  coef = lm.fit(X,Y)$coefficients
  EX = colMeans(X)
  
  for (i in 1:p ){
    for (j in 1:nsamp){
      perm = sample(1:p,p,replace=FALSE)
      k = which(perm==i)
      if (k>1){
        W = perm[1:(k-1)]
      } else{
        W = NULL
      }
      
      oi = setdiff(1:p,c(W,i))
      psi1 = X[,c(W,i),drop=FALSE] %*% as.matrix(coef[c(W,i)]) + sum(EX[oi]*coef[oi])
      
      oi = setdiff(1:p,W)
      psi2 = X[,W,drop=FALSE] %*% as.matrix(coef[W]) + sum(EX[oi]*coef[oi])
      
      S[,i] = S[,i] + (psi1 - psi2)
      
    }
    
    S[,i] = S[,i]/nsamp
    
  }
  
  return(S)
}
