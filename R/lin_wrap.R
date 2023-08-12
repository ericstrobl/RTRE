lin_wrap <- function(x,y,z,suffStat){
  
  x1=suffStat$data[,x];
  y1=suffStat$data[,y];
  z1=suffStat$data[,z,drop=FALSE];
  
  out = lin_test(x1,y1,z1)
  
  return(out$p)
  
}

lin_test <- function(x,y,z=NULL, approx="lpd4"){
  
  require(momentchi2)
  require(earth)
  
  x = normalizeData(x)
  y = normalizeData(y)
  
  r = nrow(x)
  
  xn = normalizeData(x)
  yn = normalizeData(y)
  
  
  if (length(z) > 0){
    z = normalizeData(z)
    
    res_x = as.matrix(lm.fit(z,xn)$residuals)
    res_y = as.matrix(lm.fit(z,yn)$residuals)
  } else{
    res_x = xn
    res_y = yn
  }
  
  Cxy_z = cov(res_x, res_y)
  
  
  if (ncol(xn)==1){
    approx="hbe"
  }
  
  if (approx == "perm"){
    
    Cxy_z = cov(res_x, res_y);
    Sta = r*sum(Cxy_z^2);
    
    nperm =1000;
    
    Stas = c();
    for (ps in 1:nperm){
      perm = sample(1:r,r);
      Sta_p = Sta_perm(res_x[perm,],res_y,r)
      Stas = c(Stas, Sta_p);
      
    }
    
    p = 1-(sum(Sta >= Stas)/length(Stas));
    
  } else {
    
    # Cxy_z=Cxy-Cxz%*%i_Czz%*%Czy; #less accurate for permutation testing
    Sta = r*sum(Cxy_z^2);
    
    d =expand.grid(1:ncol(xn),1:ncol(yn));
    res = res_x[,d[,1]]*res_y[,d[,2]];
    Cov = 1/r * (t(res)%*%res);
    
    if (approx == "chi2"){
      i_Cov = ginv(Cov)
      
      Sta = r * (c(Cxy_z)%*%  i_Cov %*% c(Cxy_z) );
      p = 1-pchisq(Sta, length(c(Cxy_z)));
    } else{
      
      eig_d = eigen(Cov,symmetric=TRUE);
      eig_d$values=eig_d$values[eig_d$values>0];
      
      if (approx == "gamma"){
        p=1-sw(eig_d$values,Sta);
        
      } else if (approx == "hbe") {
        
        p=1-hbe(eig_d$values,Sta);
        
      } else if (approx == "lpd4"){
        eig_d_values=eig_d$values;
        p=try(1-lpb4(eig_d_values,Sta),silent=TRUE);
        if (!is.numeric(p) | is.nan(p)){
          p=1-hbe(eig_d$values,Sta);
        }
      }
    }
  }
  
  if (p<0) p=0;
  
  out=list(p.value=p, statistic=Sta);
  return(out)
  
}