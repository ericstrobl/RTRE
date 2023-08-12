normal_reg <- function(data){
  
  p = ncol(data)-1
  S = matrix(0,nrow(data),p)
  
  coef = lm.fit(data[,1:p],data[,p+1])$coefficients
  
  for (i in 1:p){
    
    S[,i] = coef[i]*data[,i]
    
  }
  
  return(S)
  
}

