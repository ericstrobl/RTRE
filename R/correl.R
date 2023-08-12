correl <- function(data){
  
  p = ncol(data)-1
  S = matrix(0,nrow(data),p)
  
  for (i in 1:p){
    
    S[,i] = lm.fit(data[,i,drop=FALSE],data[,p+1])$coefficients[1]*data[,i]
    
  }
  
  return(S)
  
}
