#' thresholding 
#' 
#' A function
#' @param temp
#' @param tau
#' @param p
#' @return Shat
#' @export 
#' thresholding
#' 
thresholding =  function(temp,tau,p)
{
  ind = which(temp$Dval >  tau)
  
  Shat =  c()
  
  if(length(ind)==0)
  {
    return(NULL)
  }
  
  for( j in 1:length(ind))
  {
    if(p[ind[j]]==0)
    {
      Shat = c(Shat,temp$S[ind[j]])
    }
    if(p[ind[j]] > 0  && min(abs(Shat - temp$S[p[ind[j]]] ))==0   )
    {
      Shat = c(Shat,temp$S[ind[j]])
    }
  }
  
  return(Shat)
}