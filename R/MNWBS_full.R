#' MNWBS_full function
#'
#' A function to compute change points based on the  MNP method
#' @param y : a matrix containing the data. Row correspond to different time points and columns to different variables
#' @param z : a copy of the matrix y, it can be y itself
#' @param alpha : left end points of the random intervals
#' @param beta : right end points of the random intervals
#' @param h : bandwith parameter
#' @param len_tau : number of tau parameters, default is 30
#' @return est : set of estimated change points
#' @export
#' MNWBS_full
#' 
MNWBS_full =  function(y,z,alpha,beta,h,len_tau = 30)
{
  p =  dim(y)[2]
  T = dim(y)[1]
  s = 1
  e = T
  
  
  #grid = 
  #  temp1 = new_WBS(z, gam,1,T,0,NULL,NULL,1,alpha,beta,N)
  temp1 = new_MWBS(y, z,s,e,flag=0,S=NULL,Dval=NULL,pos=1,alpha,beta,8*h)
  
  
  if(length(temp1$S)>1)
  {
    Dval = temp1$Dval
    p1 =  parent(temp1)
    aux = sort(Dval,decreasing = TRUE)
    tau_grid = rev(aux[1:min(len_tau,length(Dval))])- 10^{-30}
    tau_grid =  tau_grid[which(is.na(tau_grid)==FALSE)] 
    
    
    S =  c()
    for( j in 1:length(tau_grid))
    {
      aux = thresholding(temp1,tau_grid[j],p1)
      
      if(length(aux) == 0)
        break;
      
      S[[j]] = sort(aux)
    }
  }
  if(length(temp1$S)==1)
  {
    S = list()
    S[[1]] = sort(temp1$S)
  }
  
  T= dim(y)[1]
  S = unique(S)
  if(length(S)==0)
  {
    return(NULL)
  }  
  best_count = rep(0,length(S)+1)
  
  #for(ind_p in 1:21)
  #{
  indicator = 0
  lamb = 1.8^2
  #log(sum(T))/1.5#2.5#1.5#2#2.555#
  
  h_grid = 1;# 2^seq(-3,3,length = 20)
  #c(h/3,h/2,h/1.5,h,h/1.5,2*h )
  
  #S = rev(S)
  
  
  v_n =  100
  best_count =   matrix(0, length(S)+1,v_n )
  
  #  for(ind_v in 1:v_n  )
  # {
  #  vec = runif(p)
  #  vec =  vec/sqrt(sum(v))
  min_pval = 10
  for(j in 1:length(S))
  {
    if(j == 1 && length(S[[1]])==0 )
    {
      return(B)
    }
    
    B2  =  S[[j]]
    if(j < length(S))
    {
      B1 = S[[j+1]]
    }
    if(j  == length(S))
    {
      B1 = NULL
    }
    temp = setdiff(B2,B1)
    
    st =  10^{-10}
    for(l in 1:length(temp))
    {
      eta =  temp[l]
      
      if( length(B1)==0)
      {
        eta1 = 1
        eta2 = T
      }
      if( length(B1)>0)
      {
        for(k in 1:length(B1))
        {
          if(B1[k]> eta  )
            break;
        }
        if(B1[k]> eta )
        {
          eta2 = B1[k]
          
          if(k ==1)
            eta1 = 1
          
          if(k > 1)
            eta1 = B1[k-1]+1
        }
        if(B1[k]< eta )
        {
          eta1 = B1[k]+1
          eta2 = T
        }######
      }######## if length(B1) >1
      pval  =  rep(0,length(v_n))
      vec_s =  matrix(0,v_n,p)
      for(ind_v in 1:v_n)
      {
        vec = runif(p)
        vec =  vec/sqrt(sum(vec^2))
        vec_s[ind_v,] = vec 
        if(ind_v< p+1)
        {
          vec= rep(0,p)
          vec[ind_v] =1
        }
        aux = Delta_se_t_1d(y%*%vec,eta1+1,eta2,eta,rep(1,T))  
        #val_aux[ind_v] = aux
        pval[ind_v] = exp(-2*aux^2 )
    
        
        st_aux = aux
        
      }## for ind_v
      pval_adj =  p.adjust(pval, method ="fdr" )
      #min(pval_adj)
      
      if(min_pval >  min(pval_adj))
      {
        min_pval = min(pval_adj)
      }
    }### for l
 
    if(min_pval < 0.0005)#0.001
    {
      break;
    }
    #           break;
  }##### j 
  print(min_pval)
  if(min_pval>.1)# (st <1.9)
  {
    est = NULL
    return(est)
  }
  est = S[[j]]
  
  return(est)
}
#  

