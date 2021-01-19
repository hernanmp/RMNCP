#' Delta_se_t function
#'
#' Subroutine to compute change points based on the  MNP method
#' @param y a matrix containing the data. Row correspond to different time points and columns to different variables
#' @param z a copy of the matrix y, it can be y itself
#' @param s
#' @param e 
#' @param t
#' @param h  bandwith parameter
#' @return val
#' @export
#' Delta_se_t
#' 
Delta_se_t = function(y,z,s,e,t,h)
{
  
  p =  dim(y)[2]
  # =  dim(y)[]
  
  n_st = (t-s)
  n_se = (e-s)
  n_te = (e-t)
  
  dat  = y[(s+1):t,]
  

  temp = kde(dat, gridsize = 30, eval.points = z, H = h*diag(p))
  dat  = y[(t+1):e,]
  temp2 = kde(dat, gridsize = 30, eval.points = z, H = h*diag(p))
  
  val = sqrt((n_st*n_te)/n_se)*max(abs(temp$estimate -temp2$estimate))
  
  return(val)
}