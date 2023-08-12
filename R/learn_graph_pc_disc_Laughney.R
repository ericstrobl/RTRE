learn_graph_pc_disc_Laughney <- function(data,graphT){
  
  # run PC
  p = ncol(data)
  suffStat = list()
  suffStat$data = data
  # suffStat =list(C=cor(data),n=nrow(data))
  graph = pc(suffStat,earth_wrap,0.01,p=p,u2pd="rand")@graph
  graph = as(graph,"matrix")
  
  iY = which(graph[7,]==1)
  graph[7,iY] = 0; graph[iY,7] = 1
  
  iO = which(graph!=0,arr.ind=TRUE)
  for (i in seq_len(nrow(iO))){
    if ( (graph[iO[i,1],iO[i,2]]==1) & (graph[iO[i,2],iO[i,1]]==1) ){
      if ( isAnc(graphT,iO[i,1],iO[i,2]) ){
        graph[iO[i,1],iO[i,2]]=1; graph[iO[i,2],iO[i,1]] = 0
        if (isAnc(graph,iO[i,2],iO[i,1])){ # if cyclic
          graph[iO[i,1],iO[i,2]]=0 # then remove the edge
        }
      } else{
        graph[iO[i,1],iO[i,2]]=0; graph[iO[i,2],iO[i,1]] = 1 
        if (isAnc(graph,iO[i,1],iO[i,2])){ # if cyclic
          graph[iO[i,2],iO[i,1]]=0 # then remove the edge
        }
      }
    }
  }
  
  
  return(graph)
}
