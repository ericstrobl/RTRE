RCAOn2 <- function(data,graph,Yi,nsamp=1000){
  # E includes error of Y
  
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
  
  DAG = list()
  DAG$weights = direct_effects
  DAG$graph = graph
  DAG$Y = Yi
  
  anc = isAncAll(graph,1:p,Yi)
  
  for (i in anc ){
    for (j in 1:nsamp){
      
      perm = sample(1:p,p,replace=FALSE)
      k = which(perm==i)
      if (k>1){
        W = perm[1:(k-1)]
      } else{
        W = NULL
      }
      
      psi1 = apply_RCAO(E,c(W,i),DAG,data[,Yi])
      if (k == 1){
        psi2 = mean(Y)
      } else{
        psi2 = apply_RCAO(E,W,DAG,data[,Yi])
      }
      
      gamma = psi1 - psi2
      
      S[,i] = S[,i] + gamma
      
      
    }
    
    S[,i] = S[,i]/nsamp
    
  }
  
  return(S)
}

apply_RCAO <- function(Err,W,DAG,Y){
  
  r = nrow(Err)
  boot = sample(1:r,r,replace=TRUE)
  
  Errn = Err; 
  Errn[,W] = Err[boot,W] # use same V for sample j
  Yn = sample_DAG_err2(Errn,DAG)$data[,DAG$Y] # use same V for sample j
  
  pY = -log(pmax(1-ecdf(Yn)(Y),1e-10))
  
  return(pY)
  
}