load("Laughney.RData")
data = cbind(dataf1,Y)

data = as.matrix(normalizeData(data))

Gt = matrix(0,7,7)
Gt[1,4]=1; Gt[2,3]=1; Gt[3,4]=1;
Gt[4,5]=1; Gt[1,6]=1; Gt[,7]=1; Gt[7,7] = 0

G = list()
G$graph = Gt
G$weights = Gt
err = data
for (c in 1:ncol(data)){
  pa = which(Gt[,c]>0)
  if (length(pa)>0){
    mod = lm.fit(as.matrix(cbind(data[,pa],1)),data[,c])
    
    G$weights[pa,c] = mod$coefficients[1:length(pa)]
    err[,c] = mod$residuals # regress out intervention nodes
  }
}

reps = 100
Gs = vector("list",reps)
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


for (i in 1:100){
  print(i)
  
  # bootstrap draw error terms
  ib = sample(1:nrow(data),nrow(data),replace=TRUE)
  datan = data[ib,]
  
  Yi = 7
  
  graphn = G$graph
  weightsn = G$weights
  
  #ground truth
  truth = ISVn_truth(datan,weightsn,Yi)
  GT = truth$RCEs
  cl_ans = TEs = truth$TEs
  
  ## learn graph
  DAGn = learn_graph_pc_disc_Laughney(datan,graphn)
  
  
  ## run algorithms
  print("A")
  ptm <- proc.time()
  shaps = do_reg(datan,DAGn,Yi) ##
  do_reg_res[[i]]$time = (proc.time() - ptm)[3]
  do_reg_res[[i]]$RMSE_err = sqrt(mean(  (shaps - GT)^2   ))
  do_reg_res[[i]]$RMSE_TE = sqrt(mean(  (shaps - TEs)^2   ))
  do_reg_res[[i]]$SS = root_SS(shaps,cl_ans)
  
  
  print("B")
  ptm <- proc.time()
  shaps = marg_SVn2(datan,DAGn,Yi,nsamp=200) ##
  marginal_SV_res[[i]]$time = (proc.time() - ptm)[3]
  marginal_SV_res[[i]]$RMSE_err = sqrt(mean(  (shaps - GT)^2   ))
  marginal_SV_res[[i]]$RMSE_TE = sqrt(mean(  (shaps - TEs)^2   ))
  marginal_SV_res[[i]]$SS = root_SS(shaps,cl_ans)
  
  print("C")
  ptm <- proc.time()
  shaps = conditional_SV(datan[,-Yi],datan[,Yi],nsamp=100) ##
  conditional_SV_res[[i]]$time = (proc.time() - ptm)[3]
  conditional_SV_res[[i]]$RMSE_err = sqrt(mean(  (shaps - GT)^2   ))
  conditional_SV_res[[i]]$RMSE_TE = sqrt(mean(  (shaps - TEs)^2   ))
  conditional_SV_res[[i]]$SS = root_SS(shaps,cl_ans)
  
  print("D")
  ptm <- proc.time()
  shaps = MESVn2(datan,DAGn,Yi) ##
  MESV_res[[i]]$time = (proc.time() - ptm)[3]
  MESV_res[[i]]$RMSE_err = sqrt(mean(  (shaps - GT)^2   ))
  MESV_res[[i]]$RMSE_TE = sqrt(mean(  (shaps - TEs)^2   ))
  MESV_res[[i]]$SS = root_SS(shaps,cl_ans)
  
  print("E")
  ptm <- proc.time()
  shaps = RTRE(datan,DAGn,Yi) ##
  ISV_res[[i]]$time = (proc.time() - ptm)[3]
  ISV_res[[i]]$RMSE_err = sqrt(mean(  (shaps$RCEs - GT)^2   ))
  ISV_res[[i]]$RMSE_TE = sqrt(mean(  (shaps$TEs - TEs)^2   ))
  ISV_res[[i]]$SS = root_SS(shaps$TEs,cl_ans)
  
  print("F")
  ptm <- proc.time()
  shaps = ASVn2(datan[,-Yi],datan[,Yi],DAGn,nsamp=100) ##
  ASV_res[[i]]$time = (proc.time() - ptm)[3]
  ASV_res[[i]]$RMSE_err = sqrt(mean(  (shaps - GT)^2   ))
  ASV_res[[i]]$RMSE_TE = sqrt(mean(  (shaps - TEs)^2   ))
  ASV_res[[i]]$SS = root_SS(shaps,cl_ans)
  
  print("G")
  ptm <- proc.time()
  shaps = RCAOn2(datan,DAGn,Yi,nsamp=100) ##
  RCAO_res[[i]]$time = (proc.time() - ptm)[3]
  RCAO_res[[i]]$RMSE_err = sqrt(mean(  (shaps - GT)^2   ))
  RCAO_res[[i]]$RMSE_TE = sqrt(mean(  (shaps - TEs)^2   ))
  RCAO_res[[i]]$SS = root_SS(shaps,cl_ans)
  
  print("H")
  ptm <- proc.time() 
  shaps = RCAMn2(datan,DAGn,Yi,nsamp=100) ##
  RCAM_res[[i]]$time = (proc.time() - ptm)[3]
  RCAM_res[[i]]$RMSE_err = sqrt(mean(  (shaps - GT)^2   ))
  RCAM_res[[i]]$RMSE_TE = sqrt(mean(  (shaps - TEs)^2   ))
  RCAM_res[[i]]$SS = root_SS(shaps,cl_ans)
  
  print("I")
  ptm <- proc.time()
  shaps = ICCn2(datan,DAGn,Yi,nsamp=100) ##
  ICC_res[[i]]$time = (proc.time() - ptm)[3]
  ICC_res[[i]]$RMSE_err = sqrt(mean(  (shaps - GT)^2   ))
  ICC_res[[i]]$RMSE_TE = sqrt(mean(  (shaps - TEs)^2   ))
  ICC_res[[i]]$SS = root_SS(shaps,cl_ans)
  
  print("J")
  ptm <- proc.time()
  shaps = normal_reg(datan) ##
  normal_reg_res[[i]]$time = (proc.time() - ptm)[3]
  normal_reg_res[[i]]$RMSE_err = sqrt(mean(  (shaps - GT)^2   ))
  normal_reg_res[[i]]$RMSE_TE = sqrt(mean(  (shaps - TEs)^2   ))
  normal_reg_res[[i]]$SS = root_SS(shaps,cl_ans)
  
  print("K")
  ptm <- proc.time()
  shaps = correl(datan) ##
  correl_res[[i]]$time = (proc.time() - ptm)[3]
  correl_res[[i]]$RMSE_err = sqrt(mean(  (shaps - GT)^2   ))
  correl_res[[i]]$RMSE_TE = sqrt(mean(  (shaps - TEs)^2   ))
  correl_res[[i]]$SS = root_SS(shaps,cl_ans)
  
  print("L")
  ptm <- proc.time()
  shaps = do_SVn2(datan,DAGn,Yi) ##
  do_SV_res[[i]]$time = (proc.time() - ptm)[3]
  do_SV_res[[i]]$RMSE_err = sqrt(mean(  (shaps - GT)^2   ))
  do_SV_res[[i]]$RMSE_TE = sqrt(mean(  (shaps - TEs)^2   ))
  do_SV_res[[i]]$SS = root_SS(shaps,cl_ans)
  
  save(file="Results_real_Laughney.RData",Gs,do_reg_res,marginal_SV_res,conditional_SV_res,
       MESV_res,ISV_res,ASV_res,RCAO_res,RCAM_res,ICC_res,normal_reg_res,correl_res,do_SV_res)
  
}


imax = 100
all_res = matrix(0,imax,12)


for (i in 1:imax){
  all_res[i,1] = sqrt(mean(do_reg_res[[i]]$SS))
  all_res[i,2] = sqrt(mean(marginal_SV_res[[i]]$SS))
  all_res[i,3] = sqrt(mean(conditional_SV_res[[i]]$SS))
  all_res[i,4] = sqrt(mean(MESV_res[[i]]$SS))
  all_res[i,5] = sqrt(mean(ISV_res[[i]]$SS))
  all_res[i,6] = sqrt(mean(ASV_res[[i]]$SS))
  all_res[i,7] = sqrt(mean(RCAO_res[[i]]$SS))
  all_res[i,8] = sqrt(mean(RCAM_res[[i]]$SS))
  all_res[i,9] = sqrt(mean(ICC_res[[i]]$SS))
  all_res[i,10] = sqrt(mean(normal_reg_res[[i]]$SS))
  all_res[i,11] = sqrt(mean(correl_res[[i]]$SS))
  all_res[i,12] = sqrt(mean(do_SV_res[[i]]$SS))
  
}

for (i in 1:imax){
  all_res[i,1] = sqrt(mean(do_reg_res[[i]]$RMSE_err))
  all_res[i,2] = sqrt(mean(marginal_SV_res[[i]]$RMSE_err))
  all_res[i,3] = sqrt(mean(conditional_SV_res[[i]]$RMSE_err))
  all_res[i,4] = sqrt(mean(MESV_res[[i]]$RMSE_err))
  all_res[i,5] = sqrt(mean(ISV_res[[i]]$RMSE_err))
  all_res[i,6] = sqrt(mean(ASV_res[[i]]$RMSE_err))
  all_res[i,7] = sqrt(mean(RCAO_res[[i]]$RMSE_err))
  all_res[i,8] = sqrt(mean(RCAM_res[[i]]$RMSE_err))
  all_res[i,9] = sqrt(mean(ICC_res[[i]]$RMSE_err))
  all_res[i,10] = sqrt(mean(normal_reg_res[[i]]$RMSE_err))
  all_res[i,11] = sqrt(mean(correl_res[[i]]$RMSE_err))
  all_res[i,12] = sqrt(mean(do_SV_res[[i]]$RMSE_err))
  
}

for (i in 1:imax){
  all_res[i,1] = sqrt(mean(do_reg_res[[i]]$RMSE_TE))
  all_res[i,2] = sqrt(mean(marginal_SV_res[[i]]$RMSE_TE))
  all_res[i,3] = sqrt(mean(conditional_SV_res[[i]]$RMSE_TE))
  all_res[i,4] = sqrt(mean(MESV_res[[i]]$RMSE_TE))
  all_res[i,5] = sqrt(mean(ISV_res[[i]]$RMSE_TE))
  all_res[i,6] = sqrt(mean(ASV_res[[i]]$RMSE_TE))
  all_res[i,7] = sqrt(mean(RCAO_res[[i]]$RMSE_TE))
  all_res[i,8] = sqrt(mean(RCAM_res[[i]]$RMSE_TE))
  all_res[i,9] = sqrt(mean(ICC_res[[i]]$RMSE_TE))
  all_res[i,10] = sqrt(mean(normal_reg_res[[i]]$RMSE_TE))
  all_res[i,11] = sqrt(mean(correl_res[[i]]$RMSE_TE))
  all_res[i,12] = sqrt(mean(do_SV_res[[i]]$RMSE_TE))
  
}

all_res = matrix(0,imax,10)
for (i in 1:imax){
  all_res[i,] = ISV_res[[i]]$SS
}
print(colMeans(all_res))
plot(colMeans(all_res))


for (i in 1:imax){
  all_res[i,1] = do_reg_res[[i]]$time
  all_res[i,2] = marginal_SV_res[[i]]$time
  all_res[i,3] = conditional_SV_res[[i]]$time
  all_res[i,4] = MESV_res[[i]]$time
  all_res[i,5] = ISV_res[[i]]$time
  all_res[i,6] = ASV_res[[i]]$time
  all_res[i,7] = RCAO_res[[i]]$time
  all_res[i,8] = RCAM_res[[i]]$time
  all_res[i,9] = ICC_res[[i]]$time
  all_res[i,10] = normal_reg_res[[i]]$time
  all_res[i,11] = correl_res[[i]]$time
  all_res[i,12] = do_SV_res[[i]]$time
  
}
