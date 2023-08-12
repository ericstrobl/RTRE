generate_DAG <- function(p,en){
  
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
  
  DAG$total = solve(diag(p) - DAG$weights)
  
  DAG$err = sample(1:3,p,replace=TRUE)
  
  return(DAG)
}
