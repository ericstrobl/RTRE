get_ground_truth <- function(DAG,dataT){

p = ncol(dataT$data)
ans = matrix(0,nrow(dataT$data),p-1)
totalS = solve(diag(p) - DAG$weights)
for (i in 1:(p-1)){
  
  aI = isDecAll(DAG$graph,1:(p-1),i)
  dI = isAncAll(DAG$graph,1:(p-1),DAG$Y)
  aI = setdiff(intersect(aI,dI),DAG$Y)
  for (a in aI){
    weightsT = DAG$weights
    weightsT[,a]=0
    total = solve(diag(p) - weightsT)

    if (a != i){
      #do(W) - normal
      psi1 = (dataT$E[,i] * total[i,p])
      psi2 = (dataT$E[,i] * totalS[i,p])
      # ans[,a] = ans[,a] + (dataT$E[1:10,i] * total[i,10]) - (dataT$E[1:10,i] * totalS[i,10])### first part is incorrect (dataT$E[1:10,i] * totalS[i,10]) 
    } else{
      psi1 = (0 * total[i,p])
      psi2 = (dataT$E[,i] * totalS[i,p])
      # ans[,a] = ans[,a] + (0 * total[i,10]) - (dataT$E[1:10,i] * totalS[i,10]) #psi1 - psi2
    }
    
    # if(a==3){
    #   print(c(a,i,which(DAG$graph[,i]!=0),total[c(i,which(DAG$graph[,i]!=0)),10],psi1[1]))
    # }
    
    ans[,a] = ans[,a] + (psi1-psi2)
  }
}

return(ans)

}

