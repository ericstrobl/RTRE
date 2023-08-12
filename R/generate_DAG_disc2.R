generate_DAG_disc2 <- function(p,en){
  
  N = p*p - p;
  
  DAG = list()
  DAG$Y = 0
  
  while(DAG$Y != p){
    samplesB = rbinom(N/2,1, en/(p-1) ); # sample edges
    graph = matrix(0,p,p)
    graph[upper.tri(graph, diag=FALSE)] <- samplesB; # put in edges in upper triangular
    
    ord = sample(1:p,p,replace=FALSE) # permute order
    DAG$graph = graph[ord,ord]
    
    # Ys = which( (rowSums(DAG$graph)==0) & (colSums(DAG$graph)>0) ) # no children, some observed parents
    Ys = which( colSums(DAG$graph)>0 ) # some observed parents (can have children)
    if (p %in% Ys){ # ensure that DAG$Y is the last variable
      DAG$Y = p
    }
  }
  
  weights = matrix(-0.75*runif(p^2)-0.25,p,p)
  DAG$weights = weights*DAG$graph
  DAG$weights[,DAG$Y] = DAG$weights[,DAG$Y]*sample(c(-1,1),p,replace=TRUE) # target can have negative and positive coefficients
  
  DAG$err = sample(1:3,p,replace=TRUE)
  
  discretize =  setdiff(which( (rowSums(DAG$graph) < 2)& (colSums(DAG$graph) >0)),c(DAG$Y,which(DAG$err==3))) # at least one parent, at most one child
  if (length(discretize)==1){
    DAG$discretize = discretize[sample(c(0,1),1)]
  } else if (length(discretize)>1){
    DAG$discretize = sample(discretize,sample(1:length(discretize),1))
  }
  
  graphD = DAG$graph
  weightsD = DAG$weights
  for (i in DAG$discretize){
    for (p in which(graphD[,i]!=0)){ # find parents of Xi
      ch = which(graphD[i,]!=0) # find children of Xi
      graphD[p,ch] = 1 # make children of Xi the children of Pa(Xi)
      for (c in ch){
        weightsT = DAG$weights; weightsT[,i] = 0;
        weightsD[p,c] = get_total_effects(DAG$weights,p,c) - get_total_effects(weightsT,p,c) # (total effect from p to c) - (total effect from p to c in do(Xi))
      }
    }
    graphD[i,] = 0 # no children (but keep same parents)
    weightsD[i,] = 0; weightsD[which(graphD[,i]!=0),i] = 0 # weights connected to Xi are NaN
  }
  DAG$graphO = DAG$graph; DAG$graph = graphD
  DAG$weightsO = DAG$weights; DAG$weights = weightsD
  
  # DAG$total = solve(diag(ncol(DAG$graph)) - DAG$weights)

  return(DAG)
}