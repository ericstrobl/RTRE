isAnc <- function(graph,a,b,visited=rep(FALSE,nrow(graph)))
{
  
  if (a %in% b){
    return(TRUE)
  }
  
  visited[a] = TRUE;
  
  adj = which(graph[a,] & !visited);
  
  out=FALSE;
  for (j in adj){
    out=isAnc(graph, j, b, visited);
    if(out==TRUE){
      break;
    }
  }
  
  return(out)
  
}