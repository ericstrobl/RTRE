learn_graph_pc <- function(data){
  
  # run PC
  p = ncol(data)
  suffStat = list()
  suffStat$data = data
  graph = pc(suffStat,earth_wrap,0.01,p=p,u2pd="rand")@graph
  plot(graph)
  CPDAG = as(graph,"matrix")
  res1 = my_pdag2alldags(CPDAG)
  
  n = nrow(data)
  DAGs = list()
  weights = list()
  Es = list()
  for (i in 1:length(res1)){
    DAGs[[i]] = res1[[i]]
    weights[[i]] = res1[[i]]
    E = matrix(0,n,p)
    for (j in 1:p){ #convert adjacency to weights
      pa = which(DAGs[[i]][,j]!=0)
      mod = lm.fit(cbind(data[,pa,drop=FALSE],1),data[,j])
      weights[[i]][pa,j] = mod$coefficients[1:length(pa)]
      E[,j] = mod$residuals
    }
    Es[[i]] = E
  }
  
  
  return(list(CPDAG = CPDAG, DAGs = DAGs, weights = weights, Es = Es))
}
