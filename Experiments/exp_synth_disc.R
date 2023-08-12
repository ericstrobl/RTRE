library(pcalg)

nsamps = c(1000,5000,25000)
p = 10
reps = 250

Gs = lapply(1:reps, function (.) vector("list",length(nsamps)))
do_reg_res = Gs
marginal_SV_res = Gs
conditional_SV_res = Gs
MESV_res = Gs
ISV_res = Gs
ASV_res = Gs
RCAO_res = Gs
RCAM_res = Gs
ICC_res = Gs
normal_reg_res = Gs
correl_res = Gs
do_SV_res = Gs

for (i in 1:250){
  print(i)
  DAG = generate_DAG_disc2(p,2); plot(as(DAG$graph,"graphNEL"))
  dataT = sample_DAG_start_disc2(nsamps[length(nsamps)],DAG)
  colSums(dataT$data)
  
  # ground truth treatment effects
  TEs = get_ground_truth(DAG,dataT)
  
  GT = matrix(0,nsamps[length(nsamps)],p-1)
  anc = isAncAll(DAG$graph,1:(p-1),DAG$Y)
  for (d in anc){
    GT[,d] = dataT$E[,d]*get_total_effects(DAG$weights,d,DAG$Y)
  }
  
  for (n in 1:length(nsamps)){
    
    data = dataT; 
    data$data = dataT$data[1:nsamps[n],]
    data$E = dataT$E[1:nsamps[n],]
    
    ## learn graph
    DAGn = learn_graph_pc_disc(data$data,DAG$graph)
    
    cl_ans = get_ground_truth(DAG,data)
    
    ## run algorithms
    print("A")
    ptm <- proc.time()
    shaps = do_reg(data$data,DAGn,DAG$Y) ##
    do_reg_res[[i]][[n]]$time = (proc.time() - ptm)[3]
    do_reg_res[[i]][[n]]$RMSE_err = sqrt(mean(  (shaps - GT[1:nsamps[n],])^2   ))
    do_reg_res[[i]][[n]]$RMSE_TE = sqrt(mean(  (shaps - TEs[1:nsamps[n],])^2   ))
    do_reg_res[[i]][[n]]$SS = root_SS(shaps,cl_ans)
    
    
    print("B")
    ptm <- proc.time()
    shaps = marg_SVn2(data$data,DAGn,DAG$Y,nsamp=200) ##
    marginal_SV_res[[i]][[n]]$time = (proc.time() - ptm)[3]
    marginal_SV_res[[i]][[n]]$RMSE_err = sqrt(mean(  (shaps - GT[1:nsamps[n],])^2   ))
    marginal_SV_res[[i]][[n]]$RMSE_TE = sqrt(mean(  (shaps - TEs[1:nsamps[n],])^2   ))
    marginal_SV_res[[i]][[n]]$SS = root_SS(shaps,cl_ans)
    
    print("C")
    ptm <- proc.time()
    shaps = conditional_SV(data$data[,-DAG$Y],data$data[,DAG$Y],nsamp=100) ##
    conditional_SV_res[[i]][[n]]$time = (proc.time() - ptm)[3]
    conditional_SV_res[[i]][[n]]$RMSE_err = sqrt(mean(  (shaps - GT[1:nsamps[n],])^2   ))
    conditional_SV_res[[i]][[n]]$RMSE_TE = sqrt(mean(  (shaps - TEs[1:nsamps[n],])^2   ))
    conditional_SV_res[[i]][[n]]$SS = root_SS(shaps,cl_ans)
    
    print("D")
    ptm <- proc.time()
    shaps = MESVn2(data$data,DAGn,DAG$Y) ##
    MESV_res[[i]][[n]]$time = (proc.time() - ptm)[3]
    MESV_res[[i]][[n]]$RMSE_err = sqrt(mean(  (shaps - GT[1:nsamps[n],])^2   ))
    MESV_res[[i]][[n]]$RMSE_TE = sqrt(mean(  (shaps - TEs[1:nsamps[n],])^2   ))
    MESV_res[[i]][[n]]$SS = root_SS(shaps,cl_ans)
    
    print("E")
    ptm <- proc.time()
    shaps = ISVnInt4(data$data,DAGn,DAG$Y) ##
    ISV_res[[i]][[n]]$time = (proc.time() - ptm)[3]
    ISV_res[[i]][[n]]$RMSE_err = sqrt(mean(  (shaps$RCEs - GT[1:nsamps[n],])^2   ))
    ISV_res[[i]][[n]]$RMSE_TE = sqrt(mean(  (shaps$TEs - TEs[1:nsamps[n],])^2   ))
    ISV_res[[i]][[n]]$SS = root_SS(shaps$TEs,cl_ans)
    
    print("F")
    ptm <- proc.time()
    shaps = ASVn2(data$data[,-DAG$Y],data$data[,DAG$Y],DAGn,nsamp=100) ##
    ASV_res[[i]][[n]]$time = (proc.time() - ptm)[3]
    ASV_res[[i]][[n]]$RMSE_err = sqrt(mean(  (shaps - GT[1:nsamps[n],])^2   ))
    ASV_res[[i]][[n]]$RMSE_TE = sqrt(mean(  (shaps - TEs[1:nsamps[n],])^2   ))
    ASV_res[[i]][[n]]$SS = root_SS(shaps,cl_ans)
    
    print("G")
    ptm <- proc.time()
    shaps = RCAOn2(data$data,DAGn,DAG$Y,nsamp=100) ##
    RCAO_res[[i]][[n]]$time = (proc.time() - ptm)[3]
    RCAO_res[[i]][[n]]$RMSE_err = sqrt(mean(  (shaps - GT[1:nsamps[n],])^2   ))
    RCAO_res[[i]][[n]]$RMSE_TE = sqrt(mean(  (shaps - TEs[1:nsamps[n],])^2   ))
    RCAO_res[[i]][[n]]$SS = root_SS(shaps,cl_ans)
    
    print("H")
    ptm <- proc.time() 
    shaps = RCAMn2(data$data,DAGn,DAG$Y,nsamp=100) ##
    RCAM_res[[i]][[n]]$time = (proc.time() - ptm)[3]
    RCAM_res[[i]][[n]]$RMSE_err = sqrt(mean(  (shaps - GT[1:nsamps[n],])^2   ))
    RCAM_res[[i]][[n]]$RMSE_TE = sqrt(mean(  (shaps - TEs[1:nsamps[n],])^2   ))
    RCAM_res[[i]][[n]]$SS = root_SS(shaps,cl_ans)
    
    print("I")
    ptm <- proc.time()
    shaps = ICCn2(data$data,DAGn,DAG$Y,nsamp=100) ##
    ICC_res[[i]][[n]]$time = (proc.time() - ptm)[3]
    ICC_res[[i]][[n]]$RMSE_err = sqrt(mean(  (shaps - GT[1:nsamps[n],])^2   ))
    ICC_res[[i]][[n]]$RMSE_TE = sqrt(mean(  (shaps - TEs[1:nsamps[n],])^2   ))
    ICC_res[[i]][[n]]$SS = root_SS(shaps,cl_ans)
    
    print("J")
    ptm <- proc.time()
    shaps = normal_reg(data$data) ##
    normal_reg_res[[i]][[n]]$time = (proc.time() - ptm)[3]
    normal_reg_res[[i]][[n]]$RMSE_err = sqrt(mean(  (shaps - GT[1:nsamps[n],])^2   ))
    normal_reg_res[[i]][[n]]$RMSE_TE = sqrt(mean(  (shaps - TEs[1:nsamps[n],])^2   ))
    normal_reg_res[[i]][[n]]$SS = root_SS(shaps,cl_ans)
    
    print("K")
    ptm <- proc.time()
    shaps = correl(data$data) ##
    correl_res[[i]][[n]]$time = (proc.time() - ptm)[3]
    correl_res[[i]][[n]]$RMSE_err = sqrt(mean(  (shaps - GT[1:nsamps[n],])^2   ))
    correl_res[[i]][[n]]$RMSE_TE = sqrt(mean(  (shaps - TEs[1:nsamps[n],])^2   ))
    correl_res[[i]][[n]]$SS = root_SS(shaps,cl_ans)
    
    print("L")
    ptm <- proc.time()
    shaps = do_SVn2(data$data,DAGn,DAG$Y) ##
    do_SV_res[[i]][[n]]$time = (proc.time() - ptm)[3]
    do_SV_res[[i]][[n]]$RMSE_err = sqrt(mean(  (shaps - GT[1:nsamps[n],])^2   ))
    do_SV_res[[i]][[n]]$RMSE_TE = sqrt(mean(  (shaps - TEs[1:nsamps[n],])^2   ))
    do_SV_res[[i]][[n]]$SS = root_SS(shaps,cl_ans)
    
    save(file="Results_synth.RData",Gs,do_reg_res,marginal_SV_res,conditional_SV_res,
         MESV_res,ISV_res,ASV_res,RCAO_res,RCAM_res,ICC_res,normal_reg_res,correl_res,do_SV_res)
    
  }
  
}


imax = 250
M_do_reg = matrix(0,imax,length(nsamps))
M_marginal_SV = M_do_reg
M_conditional_SV = M_do_reg
M_MESV = M_do_reg
M_ISV = M_do_reg
M_ASV = M_do_reg
M_RCAO = M_do_reg
M_RCAM = M_do_reg
M_ICC = M_do_reg
M_normal_reg = M_do_reg
M_correl = M_do_reg
M_do_SV = M_do_reg


for (i in 1:imax){
  for (n in 1:length(nsamps)){
    
    M_do_reg[i,n] = sqrt(mean(do_reg_res[[i]][[n]]$SS))
    M_marginal_SV[i,n] = sqrt(mean(marginal_SV_res[[i]][[n]]$SS))
    M_conditional_SV[i,n] = sqrt(mean(conditional_SV_res[[i]][[n]]$SS))
    M_MESV[i,n] = sqrt(mean(MESV_res[[i]][[n]]$SS))
    M_ISV[i,n] = sqrt(mean(ISV_res[[i]][[n]]$SS))
    M_ASV[i,n] = sqrt(mean(ASV_res[[i]][[n]]$SS))
    M_RCAO[i,n] = sqrt(mean(RCAO_res[[i]][[n]]$SS))
    M_RCAM[i,n] = sqrt(mean(RCAM_res[[i]][[n]]$SS))
    M_ICC[i,n] = sqrt(mean(ICC_res[[i]][[n]]$SS))
    M_normal_reg[i,n] = sqrt(mean(normal_reg_res[[i]][[n]]$SS))
    M_correl[i,n] = sqrt(mean(correl_res[[i]][[n]]$SS))
    M_do_SV[i,n] = sqrt(mean(do_SV_res[[i]][[n]]$SS))
    
  }
  
}


for (i in 1:imax){
  for (n in 1:length(nsamps)){
    
    M_do_reg[i,n] = do_reg_res[[i]][[n]]$RMSE_TE
    M_marginal_SV[i,n] = marginal_SV_res[[i]][[n]]$RMSE_TE
    M_conditional_SV[i,n] = conditional_SV_res[[i]][[n]]$RMSE_TE
    M_MESV[i,n] = MESV_res[[i]][[n]]$RMSE_TE
    M_ISV[i,n] = ISV_res[[i]][[n]]$RMSE_TE
    M_ASV[i,n] = ASV_res[[i]][[n]]$RMSE_TE
    M_RCAO[i,n] = RCAO_res[[i]][[n]]$RMSE_TE
    M_RCAM[i,n] = RCAM_res[[i]][[n]]$RMSE_TE
    M_ICC[i,n] = ICC_res[[i]][[n]]$RMSE_TE
    M_normal_reg[i,n] = normal_reg_res[[i]][[n]]$RMSE_TE
    M_correl[i,n] = correl_res[[i]][[n]]$RMSE_TE
    M_do_SV[i,n] = do_SV_res[[i]][[n]]$RMSE_TE
    
  }
  
}


for (i in 1:imax){
  for (n in 1:length(nsamps)){
    
    M_do_reg[i,n] = do_reg_res[[i]][[n]]$RMSE_err
    M_marginal_SV[i,n] = marginal_SV_res[[i]][[n]]$RMSE_err
    M_conditional_SV[i,n] = conditional_SV_res[[i]][[n]]$RMSE_err
    M_MESV[i,n] = MESV_res[[i]][[n]]$RMSE_err
    M_ISV[i,n] = ISV_res[[i]][[n]]$RMSE_err
    M_ASV[i,n] = ASV_res[[i]][[n]]$RMSE_err
    M_RCAO[i,n] = RCAO_res[[i]][[n]]$RMSE_err
    M_RCAM[i,n] = RCAM_res[[i]][[n]]$RMSE_err
    M_ICC[i,n] = ICC_res[[i]][[n]]$RMSE_err
    M_normal_reg[i,n] = normal_reg_res[[i]][[n]]$RMSE_err
    M_correl[i,n] = correl_res[[i]][[n]]$RMSE_err
    M_do_SV[i,n] = do_SV_res[[i]][[n]]$RMSE_err
    
  }
  
}


for (i in 1:imax){
  for (n in 1:length(nsamps)){
    
    M_do_reg[i,n] = do_reg_res[[i]][[n]]$time
    M_marginal_SV[i,n] = marginal_SV_res[[i]][[n]]$time
    M_conditional_SV[i,n] = conditional_SV_res[[i]][[n]]$time
    M_MESV[i,n] = MESV_res[[i]][[n]]$time
    M_ISV[i,n] = ISV_res[[i]][[n]]$time
    M_ASV[i,n] = ASV_res[[i]][[n]]$time
    M_RCAO[i,n] = RCAO_res[[i]][[n]]$time
    M_RCAM[i,n] = RCAM_res[[i]][[n]]$time
    M_ICC[i,n] = ICC_res[[i]][[n]]$time
    M_normal_reg[i,n] = normal_reg_res[[i]][[n]]$time
    M_correl[i,n] = correl_res[[i]][[n]]$time
    M_do_SV[i,n] = do_SV_res[[i]][[n]]$time
    
  }
  
}


total_ISV = matrix(0,imax,10)
total_RCAM = matrix(0,imax,10)
total_ICC = matrix(0,imax,10)
for (i in 1:imax){
  
  
  total_ISV[i,] = ISV_res[[i]][[3]]$SS
  total_ICC[i,] = ICC_res[[i]][[3]]$SS
  total_RCAM[i,] = RCAM_res[[i]][[3]]$SS
  
  
  
}
