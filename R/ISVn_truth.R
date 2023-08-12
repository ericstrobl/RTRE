ISVn_truth <- function(data,weights,Yi,nsamp=1000){
  require(earth)
  
  p = ncol(data)-1
  S = matrix(0,nrow(data),p)
  R = S
  
  # assume one graph
  graph = (weights!=0)+0
  
  # learn graph parameters
  direct_effects = weights
  
  for (tar in 1:p){
    ptar = which(graph[,tar]!=0) # find parents
    
    direct_effectsU1 = direct_effects
    direct_effectsU1[,c(tar,ptar)] = 0 # intervene on tar and ptar - so these are root vertices
    total_effectsU1 = get_total_effects(direct_effectsU1,c(tar,ptar),Yi) # can be a big inversion, just walk down instead
    
    if (length(ptar)>0){
      direct_effectsU2 = direct_effects
      direct_effectsU2[,c(ptar)] = 0 # intervene on tar and ptar - so these are root vertices
      total_effectsU2 = get_total_effects(direct_effectsU2,ptar,Yi) # can be a big inversion, just walk down instead
      psi2 = data[,c(tar,ptar),drop=FALSE] %*% as.matrix(total_effectsU1) - data[,ptar,drop=FALSE] %*% as.matrix(total_effectsU2) 
    } else{
      psi2 = data[,c(tar,ptar),drop=FALSE] %*% as.matrix(total_effectsU1) 
    }
    
    R[,tar] = psi2
    
    dI = isDecAll(graph,1:p,tar) # descendants of target
    dI = intersect(dI,isAncAll(graph,1:p,Yi)) # descendants of target and ancestors of Y
    dI = setdiff(dI,Yi) # descendants of target and ancestors of Y, not including YW
    
    for (d in dI){
      
      dataW = data
      dataW[,d] = 0 # what is d set to
      
      direct_effectsT1 = direct_effects
      direct_effectsT1[,c(d,tar,ptar)] = 0 # intervene on d, tar and ptar - so these are root vertices
      total_effects1 = get_total_effects(direct_effectsT1,c(tar,ptar),Yi) # can be a big inversion, just walk down instead
      
      if (length(ptar)>0){
        direct_effectsT2 = direct_effects
        direct_effectsT2[,c(d,ptar)] = 0 # intervene on d, tar and ptar - so these are root vertices
        total_effects2 = get_total_effects(direct_effectsT2,ptar,Yi) # can be a big inversion, just walk down instead
        psi1 = dataW[,c(tar,ptar),drop=FALSE] %*% as.matrix(total_effects1) - dataW[,ptar,drop=FALSE] %*% as.matrix(total_effects2) 
      } else{
        psi1 = dataW[,c(tar,ptar),drop=FALSE] %*% as.matrix(total_effects1) 
      }
      
      # # print(psi1[1])
      # if (d==3){
      #   print(c(d,tar,ptar,total_effects1,psi1[1]))
      # }
      # print(tar)
      S[,d] = S[,d] + (psi1 - psi2)
    }
    
  }
  
  return(list(TEs = S, RCEs = R))
  
}