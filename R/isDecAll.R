isDecAll <- function(graph,as,b){
  
  anc = c() # ancestors of b (not including b)
  
  for (a in as){
    if(isDec(graph,a,b)){
      anc = c(anc,a)
    }
  }
  
  return(anc)
}