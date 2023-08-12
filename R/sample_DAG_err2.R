sample_DAG_err2 <- function(err, DAG){
  
  G = DAG$graph
  r = nrow(G)
  
  Y = DAG$Y
  nsamps = nrow(err)
  data = err
  
  done=which(colSums(G)==0) # variables without parents
  
  stop=0;
  while (stop==0){
    for (s in done){
      ch=which(G[s,]==1) 
      for (c in ch){
        if (c %in% done){
          next
        }
        pa=which(G[,c]==1) 
        
        h=intersect(pa,done)
        if (setequal(h,pa)){ # if all parents already done
          
          A = (data[,h,drop=FALSE]%*%DAG$weights[h,c,drop=FALSE]) # do not need to include offset because the gamma error terms introduce the offset
          data[,c]= A + err[,c]
          
          done=unique(c(done, c))
        }
      }
    }
    
    if (length(done) == r){
      stop=1;
    }
  }
  
  return(list(data=data, E=err))
}