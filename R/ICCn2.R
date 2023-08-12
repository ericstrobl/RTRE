ICCn2 <- function(data,graph,Yi,nsamp=1000){
  require(earth)
  
  Y = data[,Yi]
  
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

  anc = isAncAll(graph,1:p,Yi)
  
  for (i in anc ){
    # ix = setdiff(1:p,i)
    for (j in 1:nsamp){
      
      perm = sample(1:p,p,replace=FALSE)
      k = which(perm==i)
      if (k>1){
        W = perm[1:(k-1)]
      } else{
        W = NULL
      }
      
      psi1 = lm.fit(x=E[,c(W,i),drop=FALSE],y=Y)$residuals
      psi1 = earth(x=E[,c(W,i),drop=FALSE],y=psi1^2)$fitted.values # conditional variance
      
      if (k == 1){
        psi2 = mean( (Y-mean(Y))^2 )
      } else{
        psi2 = lm.fit(x=E[,W,drop=FALSE],y=Y)$residuals
        psi2 = earth(x=E[,W,drop=FALSE],y=psi2^2)$fitted.values
      }
      
      gamma = -psi1 + psi2
      
      S[,i] = S[,i] + gamma
      
      
    }
    
    S[,i] = S[,i]/nsamp
    
  }
  
  return(S)
}