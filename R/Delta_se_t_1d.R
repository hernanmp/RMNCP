#' Delta_se_t_1d
#' 
#'  A function
#' @param y a matrix containing the data. Row correspond to different time points and columns to different variables
#' @param s
#' @param e
#' @param t
#' @param N
#' @return temp
#' @export
#' Delta_se_t_1d

Delta_se_t_1d = function(y,s,e,t,N)
{
  #T =   dim(y)[2]
  n =  dim(y)[2]
  
  n_st = sum(N[s:t])  #n*(t-s+1)
  n_se = sum(N[s:e])  #n*(e-s+1)
  n_te =sum(N[(t+1):e]) #n*(e-(t+1) +1)
  
  aux =  as.vector(y[s:t,])
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  vec_y =  as.vector(y[s:e,])
  vec_y = vec_y[which(is.na(vec_y)==FALSE)]
  Fhat_st =  temp(vec_y)# temp(grid)
  
  aux = y[(t+1):e,]
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  Fhat_te =  temp(vec_y)# temp(grid)
  
  temp =  sqrt( n_st*n_te / n_se   ) *max(abs(Fhat_te - Fhat_st  ))
  
  return(temp )
}
