sample_DAG_start_discrete <- function(nsamps, DAG){
  
  G = DAG$graph
  r = nrow(G)
  
  Y = DAG$Y
  
  done=which(colSums(G)==0) # variables without parents
  err=matrix(0,nsamps,r)
  for (i in 1:r){
    if (DAG$err[i] == 1){
      err[,i]=matrix(rt(nsamps,df=5),nsamps)
    } else if (DAG$err[i]==2){
      err[,i]=matrix(runif(nsamps,-1,1),nsamps)
    } else if (DAG$err[i]==3){
      err[,i]=matrix(rchisq(nsamps,df=3)-3,nsamps)
    }
  }
  
  data = err # create data
  
  for (c in done){ # discretize data that is already done
    if (c %in% DAG$discrete){
      # data[,c]=discretize(err[,c],DAG$cuts[c])
      data[,c]=binarize(data[,c])
    }
  }
  
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
          if (c %in% DAG$discrete){ # discretize data
            # data[,c]=discretize(data[,c],DAG$cuts[c])
            data[,c]=binarize(data[,c])
          }
          
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



discretize <- function(x,ncuts){
  require(pracma)
  percentiles = linspace(0,1,ncuts+2)
  xn = matrix(0,ncuts,length(x))
  for (c in 1:ncuts){
    xn[c,] = (x>quantile(x,percentiles[c+1]))
  }
  
  return(colSums(xn)*2-1)
  
}

binarize <- function(x){
  x[x>=0]=1
  x[x<0]=-1
  
  return(x)
  
}