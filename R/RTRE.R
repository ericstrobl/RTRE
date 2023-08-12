RTRE <- function(data,graph,Yi,do_values = rep(0,ncol(data)),exc=NULL){
  require(earth)
  
  p = ncol(data)-1
  S = matrix(0,nrow(data),p)
  R = S
  
  # learn graph parameters
  direct_effects = matrix(0,p+1,p+1)
  for (i in 1:(p+1)){
    pa = which(graph[,i]!=0)
    if (length(pa)>0){
      direct_effects[pa,i] = lm.fit(data[,pa,drop=FALSE],data[,i])$coefficients
    } 
  }
  
  anc = isAncAll(graph,1:p,Yi)
  
  for (tar in anc){
    ptar = which(graph[,tar]!=0) # find parents

    if (length(ptar)>0){
      psi2 = lm.fit(data[,c(tar,ptar),drop=FALSE],data[,Yi])$fitted.values - lm.fit(data[,ptar,drop=FALSE],data[,Yi])$fitted.values
    } else{
      psi2 = lm.fit(data[,tar,drop=FALSE],data[,Yi])$fitted.values
    } # checked

    R[,tar] = psi2
    
    dI = isDecAll(graph,1:p,tar) # descendants of target
    dI = intersect(dI,isAncAll(graph,1:p,Yi)) # descendants of target and ancestors of Y
    dI = setdiff(dI,Yi) # descendants of target and ancestors of Y, not including YW
    
    for (d in dI){
      
      if (d %in% exc){
        S[,d] = S[,d] + 0 # if d cannot be manipulated, then psi1 - psi2 = 0
      } else{
        dataW = data
        dataW[,d] = do_values[d] # what is d set to
        
        direct_effectsT1 = direct_effects
        direct_effectsT1[,c(d,tar,ptar)] = 0 # intervene on tar and ptar - so these are root vertices
        total_effects1 = get_total_effects(direct_effectsT1,c(tar,ptar),Yi) # can be a big inversion, just walk down instead
   
        if (length(ptar)>0){
          direct_effectsT2 = direct_effects
          direct_effectsT2[,c(d,ptar)] = 0 # intervene on d, tar and ptar - so these are root vertices
          total_effects2 = get_total_effects(direct_effectsT2,ptar,Yi) # can be a big inversion, just walk down instead
          psi1 = dataW[,c(tar,ptar),drop=FALSE] %*% as.matrix(total_effects1) - dataW[,ptar,drop=FALSE] %*% as.matrix(total_effects2) 
        } else{
          psi1 = dataW[,c(tar,ptar),drop=FALSE] %*% as.matrix(total_effects1) 
        }
        
        S[,d] = S[,d] + (psi1 - psi2)
      }
    }
    
  }
  
  return(list(TEs = S, RCEs = R))
  
}

train_predict <-function(xtr,ytr,xte,ic){
  xtr = as.matrix(xtr); xte = as.matrix(xte)
  
  sds = apply(xtr,2,sd)
  ir = which(sds<1E-6)
  ic = setdiff(ic,ir) # remove variables with zero variance
  
  off = rep(1,length(ytr))
  xtr = cbind(xtr[,ic,drop=FALSE],off); xte = cbind(xte[,ic,drop=FALSE],off) # only include variables with variances, add in offset
  if (ncol(xtr)>0){
    coef = as.matrix(lm.fit(x=as.matrix(xtr),y=ytr)$coefficients)
    # print(coef)
    return( xte %*% coef ) # E(theta|X,Pre(X),W)
  } else{
    return(rep(mean(ytr),length(ytr)))
  }
}
