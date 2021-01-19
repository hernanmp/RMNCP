#' new_MWBS function
#'
#' Sub-routine compute change points based on the  MNP method
#' @param y a matrix containing the data. Row correspond to different time points and columns to different variables
#' @param z a copy of the matrix y, it can be y itself
#' @param s
#' @param e
#' @param flag  
#' @param S
#' @param Dval
#' @param pos
#' @param alpha
#' @param beta
#' @param h
#' @return S
#' @return Dval
#' @return pos
#' @export
#' new_MWBS
#' 
new_MWBS =  function(y, z,s,e,flag=0,S=NULL,Dval=NULL,pos=1,alpha,beta,h)
{
  #print("pos")
  #print(pos)
  if(s> e)
  {
    S =  c(S,NULL)
    Dval =  c(Dval, NULL)
    pos  = NULL#c(pos,NULL)
    return(list(S=S,pos = pos,Dval = Dval))
  }
  
  alpha_new =  pmax(alpha,s)
  beta_new = pmin(beta,e)
  p =  dim(y)[2]
  # 
  #  ind = which( beta_new- alpha_new > 1+  3*gam)
  #  alpha_new =  alpha_new[ind]
  #  beta_new=   beta_new[ind]
  # # 
  ind = which( beta_new- alpha_new >2*max(p,h^{-p} ))#3*max(p,h^{-p} ) )
  alpha_new =  alpha_new[ind]
  beta_new=   beta_new[ind]
  M =  length(alpha_new)
  
  xi = 1/8
  alpha_new2 =alpha_new
  beta_new2  =  beta_new
  alpha_new =   ceiling((1-xi)*alpha_new2+  xi*beta_new2)
  beta_new =    ceiling((1-xi)*beta_new2 +  xi*alpha_new2)
  ind = which( beta_new- alpha_new >1)
  alpha_new =  alpha_new[ind]
  beta_new=   beta_new[ind]
  M =  length(alpha_new)
  
  # print(S)
  if(M==0 && e-s >2*max(p,h^{-p} ))#3*max(p,h^{-p} )   )
  {
    alpha_new =s
    beta_new =e
  }
  #  beta_new=   beta_new[ind]
  # # 
  #  ind = whic
  if(M ==  0 ||  pos[length(pos)]>7)
  {
    S =  c(S,NULL)
    Dval =  c(Dval, NULL)
    pos  = NULL#c(pos,NULL)
    return(list(S=S,pos = pos,Dval = Dval))
  }
  
  b  =  rep(0,M)
  a =  rep(0,M)
  
  #print(beta_new)
  #print(alpha_new)
  #  print()
  for( m in 1:M  )
  {
    temp  =  rep(0,beta_new[m]-alpha_new[m]+1)
    for(t  in    (alpha_new[m]+   max(h^{-p},p)):(beta_new[m]-max(h^{-p},p) )  )
    {
      temp[t-(alpha_new[m]) ] =  Delta_se_t(y,z,alpha_new[m],beta_new[m],t,h)
      #Delta_se_t(y,alpha_new[m],beta_new[m],t,N)
    }
    best_ind  =  which.max(temp)
    a[m] =  alpha_new[m] +  best_ind
    b[m] =   temp[best_ind]
  }
  best_ind =  which.max(b)
  
  if(abs(a[best_ind]  - s )  <= p   ||  abs(a[best_ind]  - e )  <= p  )
  {
    S =  c(S,NULL)
    Dval =  c(Dval, NULL)
    pos  = NULL#c(pos,NULL)
    return(list(S=S,pos = pos,Dval = Dval))
  }
  #  print(b)
  #  if(b[best_ind]  <tau  )
  # {
  #  return(S)
  #}
  # print(S)
  #best_t =   s+ best_t 
  pos1 =  pos
  pos2 = pos
  pos1[length(pos)] = pos[length(pos)]+1
  pos2[length(pos)] = pos[length(pos)]+1
  
  #print(c(a[best_ind],b[best_ind]))
  temp2 = new_MWBS(y, z,a[best_ind]+p,e,flag,S,Dval,pos2,alpha,beta,h)
  #new_WBS(y,gam,a[best_ind]+1,e,flag, S,Dval,pos2,alpha, beta,N)
  temp1 = new_MWBS(y, z,s,a[best_ind]-p,flag,S,Dval,pos1,alpha,beta,h)
  #new_WBS(y,gam,s,a[best_ind]-1,flag, S,Dval,pos1,alpha, beta,N)
  S1 = temp1$S 
  Dval1 = temp1$Dval     
  pos1 = temp1$pos  
  S2 = temp2$S 
  Dval2 = temp2$Dval 
  pos2 = temp2$pos 
  
  S =  c(S,a[best_ind])
  Dval = c(Dval,b[best_ind])
  
  
  S  =   c(S,S1,S2)
  Dval =  c(Dval,Dval1,Dval2)
  pos = c(pos,pos1,pos2)
  
  return(list(S=S,Dval = Dval,pos=pos))
  # S  =   c(S,S1,S2)
  # 
  # 
  # return(S)
  
}